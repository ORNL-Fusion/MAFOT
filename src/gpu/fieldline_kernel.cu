// CUDA implementation of MAFOT field-line tracing.
// Compiled by nvcc.  No blitz++ or MPI headers.
//
// Architecture:
//   CPU pre-samples getBfield() onto a uniform 3-D (R,phi,Z) grid.
//   Each GPU thread traces one field line using trilinear interpolation
//   from that grid — so all equilibrium backends (EFIT, M3DC1, GPEC, VMEC,
//   user-supplied netCDF) are supported transparently.

#include <fieldline_kernel.cuh>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

// ---------------------------------------------------------------------------
// CUDA error-checking macro
// ---------------------------------------------------------------------------
#define CUDA_CHECK(call)                                                      \
    do {                                                                       \
        cudaError_t _e = (call);                                               \
        if(_e != cudaSuccess) {                                                \
            fprintf(stderr, "CUDA error %s:%d  %s\n",                         \
                    __FILE__, __LINE__, cudaGetErrorString(_e));                \
            return -1;                                                          \
        }                                                                      \
    } while(0)

// ---------------------------------------------------------------------------
// Number of CUDA-capable devices (-1 if the CUDA runtime errors, e.g. no
// driver). Lets the callers fail clearly when -g is given without a GPU.
// ---------------------------------------------------------------------------
int gpu_device_count(void)
{
    int n = 0;
    if(cudaGetDeviceCount(&n) != cudaSuccess) return -1;
    return n;
}

// ---------------------------------------------------------------------------
// Device-side mirror of FieldGrid3D — holds device pointers.
// One of these is copied into device constant memory before kernel launch.
// ---------------------------------------------------------------------------
struct DevGrid
{
    const double* BR;
    const double* Bphi;
    const double* BZ;
    const double* psi;       // 2-D [NR*NZ]

    const double* wall_R;
    const double* wall_Z;
    int    Nwall;

    int    NR, Nphi, NZ;
    double Rmin, dR;
    double phimin, dphi;
    double Zmin,  dZ;
    double phimax;           // for phi-periodicity clamping

    double bndy[4];
    int    simpleBndy;

    double v_par, v_radial, v_tor;   // prescribed drift velocities [m/s]; v_radial=v_tor=0 -> pure field line

    int    sigma, Zq;                // relativistic guiding-center drift (sigma==0 -> off)
    double GAMMA, eps0, Ix, R0;      // drift constants of motion (from host PARTICLE)
};

__constant__ DevGrid c_grid;   // constant memory copy

// ============================================================
//  Device helper functions
// ============================================================

static __device__ __forceinline__
double clamp_d(double v, double lo, double hi)
{
    return fmax(lo, fmin(hi, v));
}

// ---------------------------------------------------------------------------
// Trilinear interpolation of one B-field component on the (R,phi,Z) grid.
// Handles phi-periodicity and clamps out-of-bounds R/Z to edge values.
// ---------------------------------------------------------------------------
static __device__
double interp3D(const double* __restrict__ F,
                double R, double phi_deg, double Z)
{
    // --- R index ---
    double fR = (R - c_grid.Rmin) / c_grid.dR;
    fR = clamp_d(fR, 0.0, (double)(c_grid.NR - 1));
    int iR = (int)fR;
    if(iR >= c_grid.NR - 1) iR = c_grid.NR - 2;
    double tR = fR - iR;

    // --- phi index (periodic) ---
    double fphi;
    if(c_grid.Nphi <= 1) {
        fphi = 0.0;
    } else {
        // normalise phi to [phimin, phimax)
        double phi_range = c_grid.phimax - c_grid.phimin;
        double pmod = fmod(phi_deg - c_grid.phimin, phi_range);
        if(pmod < 0) pmod += phi_range;
        fphi = pmod / c_grid.dphi;
        if(fphi >= c_grid.Nphi) fphi -= c_grid.Nphi;
    }
    int iphi = (int)fphi;
    // Only clamp the upper edge for true 3-D grids. For Nphi==1 (axisymmetric),
    // Nphi-2 == -1, which would index the wrong R-row (or out of bounds) and
    // corrupt every B-field lookup -> field line crosses flux surfaces / diverges.
    if(c_grid.Nphi > 1 && iphi >= c_grid.Nphi - 1) iphi = c_grid.Nphi - 2;
    double tphi = (c_grid.Nphi > 1) ? (fphi - iphi) : 0.0;

    // --- Z index ---
    double fZ = (Z - c_grid.Zmin) / c_grid.dZ;
    fZ = clamp_d(fZ, 0.0, (double)(c_grid.NZ - 1));
    int iZ = (int)fZ;
    if(iZ >= c_grid.NZ - 1) iZ = c_grid.NZ - 2;
    double tZ = fZ - iZ;

    int Nphi = c_grid.Nphi;
    int NZ   = c_grid.NZ;
    int iphi1 = (Nphi > 1) ? ((iphi + 1) % Nphi) : 0;

    // 8 corner indices
    #define IDX(r,p,z) ((long)(r)*Nphi*NZ + (long)(p)*NZ + (z))
    double c000 = F[IDX(iR,   iphi,  iZ  )];
    double c001 = F[IDX(iR,   iphi,  iZ+1)];
    double c010 = F[IDX(iR,   iphi1, iZ  )];
    double c011 = F[IDX(iR,   iphi1, iZ+1)];
    double c100 = F[IDX(iR+1, iphi,  iZ  )];
    double c101 = F[IDX(iR+1, iphi,  iZ+1)];
    double c110 = F[IDX(iR+1, iphi1, iZ  )];
    double c111 = F[IDX(iR+1, iphi1, iZ+1)];
    #undef IDX

    // Trilinear interpolation
    double v = c000*(1-tR)*(1-tphi)*(1-tZ)
             + c001*(1-tR)*(1-tphi)*   tZ
             + c010*(1-tR)*   tphi *(1-tZ)
             + c011*(1-tR)*   tphi *   tZ
             + c100*   tR *(1-tphi)*(1-tZ)
             + c101*   tR *(1-tphi)*   tZ
             + c110*   tR *   tphi *(1-tZ)
             + c111*   tR *   tphi *   tZ;
    return v;
}

// ---------------------------------------------------------------------------
// Bilinear interpolation of the 2-D psi grid
// ---------------------------------------------------------------------------
static __device__
double interp2D_psi(double R, double Z)
{
    double fR = (R - c_grid.Rmin) / c_grid.dR;
    fR = clamp_d(fR, 0.0, (double)(c_grid.NR - 1));
    int iR = (int)fR; if(iR >= c_grid.NR-1) iR = c_grid.NR-2;
    double tR = fR - iR;

    double fZ = (Z - c_grid.Zmin) / c_grid.dZ;
    fZ = clamp_d(fZ, 0.0, (double)(c_grid.NZ - 1));
    int iZ = (int)fZ; if(iZ >= c_grid.NZ-1) iZ = c_grid.NZ-2;
    double tZ = fZ - iZ;

    int NZ = c_grid.NZ;
    double v = c_grid.psi[(long)iR    *NZ + iZ  ]*(1-tR)*(1-tZ)
             + c_grid.psi[(long)iR    *NZ + iZ+1]*(1-tR)*   tZ
             + c_grid.psi[(long)(iR+1)*NZ + iZ  ]*   tR *(1-tZ)
             + c_grid.psi[(long)(iR+1)*NZ + iZ+1]*   tR *   tZ;
    return v;
}

// ---------------------------------------------------------------------------
// Winding-number point-in-polygon test (port of outofRealBndy from mafot.hxx).
// Returns true if (R,Z) is OUTSIDE the wall polygon.
// ---------------------------------------------------------------------------
static __device__
bool outsideWall(double R, double Z)
{
    int Nwall = c_grid.Nwall;
    if(Nwall < 3) return false;   // no wall loaded → never outside

    int wn = 0;
    double x2 = c_grid.wall_R[0];
    double y2 = c_grid.wall_Z[0];
    bool startOver = (y2 >= Z);

    for(int i = 1; i < Nwall; i++) {
        if(x2 == c_grid.wall_R[i] && y2 == c_grid.wall_Z[i]) continue;

        double x1 = x2, y1 = y2;
        x2 = c_grid.wall_R[i];
        y2 = c_grid.wall_Z[i];

        // Edge cases: point on edge → treat as inside
        if((y1==y2) && (Z==y1) && (((x1<=R)&&(R<=x2))||((x2<=R)&&(R<=x1)))) return false;
        if((x1==x2) && (R==x1) && (((y1<=Z)&&(Z<=y2))||((y2<=Z)&&(Z<=y1)))) return false;
        double a = (R-x1)*(y2-y1)-(Z-y1)*(x2-x1);
        if(a <= 1e-15 && a >= -1e-15) return false;

        bool endOver = (y2 >= Z);
        if(startOver != endOver) {
            if((y2-Z)*(x2-x1) <= (y2-y1)*(x2-R)) { if(endOver)  wn++; }
            else                                    { if(!endOver) wn--; }
        }
        startOver = endOver;
    }
    return (wn == 0);   // wn==0 → outside
}

// ---------------------------------------------------------------------------
// Boundary check: returns true if the field line should terminate here.
// ---------------------------------------------------------------------------
static __device__ __forceinline__
bool outofBndy(double R, double Z)
{
    if(isnan(R) || isnan(Z) || isinf(R) || isinf(Z)) return true;
    if(c_grid.simpleBndy == 1) {
        return (R < c_grid.bndy[0] || R > c_grid.bndy[1] ||
                Z < c_grid.bndy[2] || Z > c_grid.bndy[3]);
    }
    return outsideWall(R, Z);
}

// ---------------------------------------------------------------------------
// Field-line ODE right-hand side: dR/dphi, dZ/dphi
// phi_deg is current phi in degrees (right-handed).
// Returns false if Bphi is too small (degenerate field).
// ---------------------------------------------------------------------------
static __device__ __forceinline__
bool dgls(double R, double Z, double phi_deg,
          double& dRdphi, double& dZdphi)
{
    double BR   = interp3D(c_grid.BR,   R, phi_deg, Z);
    double BZ   = interp3D(c_grid.BZ,   R, phi_deg, Z);
    double Bphi = interp3D(c_grid.Bphi, R, phi_deg, Z);

    if(fabs(Bphi) < 1e-30) return false;
    // The field-line ODE (dR/dphi = R*BR/Bphi) is per-RADIAN, but this integrator advances phi in
    // DEGREES (interp3D/output/boundary all use degrees), so scale by pi/180 for the per-degree
    // slope (the CPU does this by converting dpinit[deg]->rad, particle_class.hxx).
    const double DEG2RAD = 0.017453292519943295;
    if(c_grid.v_radial != 0.0 || c_grid.v_tor != 0.0)
    {
        // Prescribed drift-velocity field v = v_par*b_hat + v_radial*r_hat + v_tor*phi_hat,
        // reparametrized by phi. r_hat = (-BZ,0,BR)/Bp (poloidal flux-surface normal, matches
        // HEAT fluxSurfNorms). v_radial=v_tor=0 reduces exactly to the field-line slope below.
        double Bp   = sqrt(BR*BR + BZ*BZ);
        double modB = sqrt(BR*BR + BZ*BZ + Bphi*Bphi);
        double vR   = c_grid.v_par*BR/modB - c_grid.v_radial*BZ/Bp;
        double vZ   = c_grid.v_par*BZ/modB + c_grid.v_radial*BR/Bp;
        double vphi = c_grid.v_par*Bphi/modB + c_grid.v_tor;
        dRdphi = R * vR / vphi * DEG2RAD;
        dZdphi = R * vZ / vphi * DEG2RAD;
    }
    else
    {
        dRdphi = R * BR / Bphi * DEG2RAD;
        dZdphi = R * BZ / Bphi * DEG2RAD;
    }

    // Relativistic guiding-center drift (grad-B + curvature), added to dZ/dphi.  Mirrors the
    // sigma!=0 branch of PARTICLE::dgls with no E-field / no sheath (gamma == GAMMA).  GAMMA, eps0,
    // Ix are constants of motion computed on the host PARTICLE.
    if(c_grid.sigma != 0)
    {
        double Sarg = c_grid.eps0*(c_grid.GAMMA*c_grid.GAMMA - 1.0) - 2.0*c_grid.R0*c_grid.Ix/R;
        if(Sarg < 0.0) return false;   // canonical-momentum constraint violated -> terminate line
        double S = sqrt(Sarg);
        dZdphi += -((double)c_grid.sigma/(double)c_grid.Zq) * (R*S + c_grid.R0*c_grid.Ix/S) * DEG2RAD;
    }
    return true;
}

// ---------------------------------------------------------------------------
// One 4th-order Runge-Kutta step.
// (R, Z, phi_deg) updated in-place; phi advanced by dphi_deg.
// Returns false on field degeneracy.
// ---------------------------------------------------------------------------
static __device__
bool rk4_step(double& R, double& Z, double& phi_deg, double dphi_deg)
{
    double k1R, k1Z, k2R, k2Z, k3R, k3Z, k4R, k4Z;
    double h2 = dphi_deg * 0.5;

    if(!dgls(R,             Z,             phi_deg,      k1R, k1Z)) return false;
    if(!dgls(R + h2*k1R,   Z + h2*k1Z,   phi_deg + h2, k2R, k2Z)) return false;
    if(!dgls(R + h2*k2R,   Z + h2*k2Z,   phi_deg + h2, k3R, k3Z)) return false;
    if(!dgls(R +    dphi_deg*k3R, Z +    dphi_deg*k3Z, phi_deg+dphi_deg, k4R, k4Z)) return false;

    R       += (dphi_deg / 6.0) * (k1R + 2*k2R + 2*k3R + k4R);
    Z       += (dphi_deg / 6.0) * (k1Z + 2*k2Z + 2*k3Z + k4Z);
    phi_deg += dphi_deg;
    return true;
}

// ============================================================
//  Global kernels
// ============================================================

// ---------------------------------------------------------------------------
// Laminar kernel — one thread per point
// ---------------------------------------------------------------------------
__global__
void k_laminar(const double* __restrict__ init_R,
               const double* __restrict__ init_Z,
               const double* __restrict__ init_phi,
               double* out_ntor,
               double* out_Lc,
               double* out_psimin,
               double* out_psiav,
               int N,
               double dpinit_deg,
               int    itt,
               int    MapDirection)
{
    int i = (int)(blockIdx.x * blockDim.x) + threadIdx.x;
    if(i >= N) return;

    double Lc_total    = 0.0;
    double psimin_tot  = 10.0;
    double psiav_tot   = 0.0;
    double phi_total   = 0.0;   // |phi| accumulated across both directions
    int    n_psi       = 0;

    // Trace in each direction
    for(int dir = -1; dir <= 1; dir += 2)
    {
        if(dir == -1 && MapDirection > 0) continue;
        if(dir ==  1 && MapDirection < 0) continue;

        double R       = init_R  [i];
        double Z       = init_Z  [i];
        double phi     = init_phi[i];
        double step    = dir * dpinit_deg;
        double phi0    = phi;
        double Lc      = 0.0;
        double psimin  = 10.0;

        // max steps = itt full toroidal turns
        int max_steps = (int)((double)itt * 360.0 / fabs(step) + 0.5);

        for(int k = 0; k < max_steps; k++)
        {
            double R_old = R, Z_old = Z, phi_old = phi;

            if(!rk4_step(R, Z, phi, step)) break;
            if(outofBndy(R, Z)) break;

            // Arc length (only inside real wall when simpleBndy==0)
            bool inside_wall = (c_grid.simpleBndy == 1) ? true : !outsideWall(R, Z);
            if(inside_wall) {
                double dR = R - R_old, dZ = Z - Z_old, dphi_rad = (phi - phi_old)*0.017453292519943;
                double Ravg = 0.5*(R + R_old);
                Lc += sqrt(dR*dR + dZ*dZ + Ravg*Ravg*dphi_rad*dphi_rad);
            }

            // Psi tracking
            double psi_v = interp2D_psi(R, Z);
            if(psi_v < psimin) psimin = psi_v;
            psiav_tot += psi_v;
            n_psi++;
        }

        Lc_total   += Lc;
        phi_total  += fabs(phi - phi0);
        if(psimin < psimin_tot) psimin_tot = psimin;
    }

    out_ntor  [i] = phi_total / 360.0;
    out_Lc    [i] = Lc_total;
    out_psimin[i] = psimin_tot;
    out_psiav [i] = (n_psi > 0) ? (psiav_tot / n_psi) : 0.0;
}

// ---------------------------------------------------------------------------
// Structure kernel — one thread per field line, one direction per launch
// direction: +1 or -1
// ---------------------------------------------------------------------------
__global__
void k_structure(const double* __restrict__ init_R,
                 const double* __restrict__ init_Z,
                 const double* __restrict__ init_phi,
                 StructureStep* results,    // [N * max_steps]
                 int*           nsteps_out, // [N]
                 int N, int max_steps,
                 double dpinit_deg, int nstep_out, int itt,
                 int direction)             // +1 or -1
{
    int i = (int)(blockIdx.x * blockDim.x) + threadIdx.x;
    if(i >= N) return;

    double R    = init_R  [i];
    double Z    = init_Z  [i];
    double phi  = init_phi[i];
    double step = direction * dpinit_deg;

    // Pre-mark all output steps invalid
    for(int j = 0; j < max_steps; j++)
        results[(long)i*max_steps + j].valid = 0;

    double Lc     = 0.0;
    int    nout   = 0;

    for(int j = 0; j < itt && nout < max_steps; j++)
    {
        // Advance nstep_out integration sub-steps (each dpinit_deg)
        bool   hit         = false;
        double lcstep_last = 0.0;        // arc length of the final sub-step
        double R_pre = R, Z_pre = Z;     // position just before the final sub-step
        for(int k = 0; k < nstep_out; k++)
        {
            double R_old = R, Z_old = Z, phi_old = phi;
            if(!rk4_step(R, Z, phi, step)) { hit = true; break; }
            if(outofBndy(R, Z))            { hit = true; break; }

            // Arc length of this sub-step
            double dR = R - R_old, dZ = Z - Z_old;
            double dphi_rad = (phi - phi_old)*0.017453292519943;
            double Ravg = 0.5*(R + R_old);
            double lcstep = sqrt(dR*dR + dZ*dZ + Ravg*Ravg*dphi_rad*dphi_rad);

            bool inside_wall = (c_grid.simpleBndy == 1) ? true : !outsideWall(R, Z);
            if(inside_wall) Lc += lcstep;

            lcstep_last = lcstep;        // remember the last sub-step
            R_pre = R_old; Z_pre = Z_old;
        }

        // Record output step
        double psi_v = interp2D_psi(R, Z);
        // dpsidLc = local gradient over the final sub-step, matching the CPU
        // (particle_class.hxx: dpsidLc = sign(dx)*(psi-psiold)/lcstep) — NOT
        // divided by the cumulative connection length.
        double psi_pre   = interp2D_psi(R_pre, Z_pre);
        double dpsidLc_v = (lcstep_last > 1e-30)
                         ? (double)direction * (psi_v - psi_pre) / lcstep_last
                         : 0.0;

        StructureStep& s = results[(long)i*max_steps + nout];
        s.R       = R;
        s.Z       = Z;
        s.phi     = phi;
        s.psi     = psi_v;
        s.Lc      = Lc;
        s.dpsidLc = dpsidLc_v;
        s.valid   = 1;
        nout++;

        if(hit) break;
    }

    nsteps_out[i] = nout;
}

// ============================================================
//  Host wrapper functions
// ============================================================

// Copy the FieldGrid3D data to device and build a DevGrid in constant memory.
// Returns 0 on success.
static int upload_grid(const FieldGrid3D* h, const GPUTraceParams& params,
                       double** d_BR, double** d_Bphi, double** d_BZ,
                       double** d_psi,
                       double** d_wall_R, double** d_wall_Z)
{
    long Nb  = (long)h->NR * h->Nphi * h->NZ * sizeof(double);
    long Np  = (long)h->NR * h->NZ           * sizeof(double);
    long Nw  = (long)h->Nwall                * sizeof(double);

    CUDA_CHECK(cudaMalloc(d_BR,   Nb)); CUDA_CHECK(cudaMemcpy(*d_BR,   h->BR,   Nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(d_Bphi, Nb)); CUDA_CHECK(cudaMemcpy(*d_Bphi, h->Bphi, Nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(d_BZ,   Nb)); CUDA_CHECK(cudaMemcpy(*d_BZ,   h->BZ,   Nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(d_psi,  Np)); CUDA_CHECK(cudaMemcpy(*d_psi,  h->psi_norm, Np, cudaMemcpyHostToDevice));

    if(h->Nwall > 0 && h->wall_R && h->wall_Z) {
        CUDA_CHECK(cudaMalloc(d_wall_R, Nw)); CUDA_CHECK(cudaMemcpy(*d_wall_R, h->wall_R, Nw, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMalloc(d_wall_Z, Nw)); CUDA_CHECK(cudaMemcpy(*d_wall_Z, h->wall_Z, Nw, cudaMemcpyHostToDevice));
    } else {
        *d_wall_R = nullptr;
        *d_wall_Z = nullptr;
    }

    DevGrid dg;
    dg.BR    = *d_BR;   dg.Bphi = *d_Bphi; dg.BZ   = *d_BZ;
    dg.psi   = *d_psi;
    dg.wall_R = *d_wall_R; dg.wall_Z = *d_wall_Z; dg.Nwall = h->Nwall;
    dg.NR = h->NR; dg.Nphi = h->Nphi; dg.NZ = h->NZ;
    dg.Rmin   = h->Rmin;   dg.dR   = h->dR;
    dg.phimin = h->phimin; dg.dphi = h->dphi; dg.phimax = h->phimax;
    dg.Zmin   = h->Zmin;   dg.dZ   = h->dZ;
    for(int k = 0; k < 4; k++) dg.bndy[k] = h->bndy[k];
    dg.simpleBndy = h->simpleBndy;
    dg.v_par = params.v_par; dg.v_radial = params.v_radial; dg.v_tor = params.v_tor;
    dg.sigma = params.sigma; dg.Zq = params.Zq;
    dg.GAMMA = params.GAMMA; dg.eps0 = params.eps0; dg.Ix = params.Ix; dg.R0 = params.R0;

    CUDA_CHECK(cudaMemcpyToSymbol(c_grid, &dg, sizeof(DevGrid)));
    return 0;
}

// ---------------------------------------------------------------------------
// gpu_trace_laminar
// ---------------------------------------------------------------------------
int gpu_trace_laminar(const FieldlineInit* init,
                      LaminarResult*       results,
                      int                  N,
                      const FieldGrid3D*   grid,
                      const GPUTraceParams& params)
{
    // --- Upload grid ---
    double *d_BR, *d_Bphi, *d_BZ, *d_psi, *d_wall_R, *d_wall_Z;
    if(upload_grid(grid, params, &d_BR, &d_Bphi, &d_BZ, &d_psi, &d_wall_R, &d_wall_Z) != 0)
        return -1;

    // --- Upload initial conditions ---
    double *d_initR, *d_initZ, *d_initPhi;
    double *h_R   = new double[N];
    double *h_Z   = new double[N];
    double *h_phi = new double[N];
    for(int i = 0; i < N; i++) { h_R[i]=init[i].R; h_Z[i]=init[i].Z; h_phi[i]=init[i].phi; }

    CUDA_CHECK(cudaMalloc(&d_initR,   N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_initZ,   N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_initPhi, N*sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_initR,   h_R,   N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_initZ,   h_Z,   N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_initPhi, h_phi, N*sizeof(double), cudaMemcpyHostToDevice));
    delete[] h_R; delete[] h_Z; delete[] h_phi;

    // --- Allocate output on device ---
    double *d_ntor, *d_Lc, *d_psimin, *d_psiav;
    CUDA_CHECK(cudaMalloc(&d_ntor,   N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Lc,     N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_psimin, N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_psiav,  N*sizeof(double)));

    // --- Launch ---
    int threads = 256;
    int blocks  = (N + threads - 1) / threads;

    printf("GPU: launching laminar kernel (%d threads, %d blocks) ...\n", threads, blocks);
    fflush(stdout);

    k_laminar<<<blocks, threads>>>(d_initR, d_initZ, d_initPhi,
                                    d_ntor, d_Lc, d_psimin, d_psiav,
                                    N, params.dpinit, params.itt, params.MapDirection);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    printf("GPU: kernel complete, retrieving results ...\n");
    fflush(stdout);

    // --- Retrieve ---
    double *h_ntor   = new double[N];
    double *h_Lc     = new double[N];
    double *h_psimin = new double[N];
    double *h_psiav  = new double[N];

    CUDA_CHECK(cudaMemcpy(h_ntor,   d_ntor,   N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_Lc,     d_Lc,     N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_psimin, d_psimin, N*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_psiav,  d_psiav,  N*sizeof(double), cudaMemcpyDeviceToHost));

    for(int i = 0; i < N; i++) {
        results[i].ntor   = h_ntor  [i];
        results[i].Lc     = h_Lc    [i];
        results[i].psimin = h_psimin[i];
        results[i].psiav  = h_psiav [i];
    }

    delete[] h_ntor; delete[] h_Lc; delete[] h_psimin; delete[] h_psiav;

    // --- Free device memory ---
    cudaFree(d_BR); cudaFree(d_Bphi); cudaFree(d_BZ); cudaFree(d_psi);
    if(d_wall_R) cudaFree(d_wall_R);
    if(d_wall_Z) cudaFree(d_wall_Z);
    cudaFree(d_initR); cudaFree(d_initZ); cudaFree(d_initPhi);
    cudaFree(d_ntor);  cudaFree(d_Lc); cudaFree(d_psimin); cudaFree(d_psiav);

    return 0;
}

// ---------------------------------------------------------------------------
// gpu_trace_structure
// ---------------------------------------------------------------------------
int gpu_trace_structure(const FieldlineInit* init,
                        int                  N,
                        int                  max_steps,
                        StructureStep*       results_bwd,
                        StructureStep*       results_fwd,
                        int*                 steps_bwd,
                        int*                 steps_fwd,
                        const FieldGrid3D*   grid,
                        const GPUTraceParams& params)
{
    // --- Upload grid ---
    double *d_BR, *d_Bphi, *d_BZ, *d_psi, *d_wall_R, *d_wall_Z;
    if(upload_grid(grid, params, &d_BR, &d_Bphi, &d_BZ, &d_psi, &d_wall_R, &d_wall_Z) != 0)
        return -1;

    // --- Upload initial conditions ---
    double *d_initR, *d_initZ, *d_initPhi;
    double *h_R   = new double[N];
    double *h_Z   = new double[N];
    double *h_phi = new double[N];
    for(int i = 0; i < N; i++) { h_R[i]=init[i].R; h_Z[i]=init[i].Z; h_phi[i]=init[i].phi; }
    CUDA_CHECK(cudaMalloc(&d_initR,   N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_initZ,   N*sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_initPhi, N*sizeof(double)));
    CUDA_CHECK(cudaMemcpy(d_initR,   h_R,   N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_initZ,   h_Z,   N*sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_initPhi, h_phi, N*sizeof(double), cudaMemcpyHostToDevice));
    delete[] h_R; delete[] h_Z; delete[] h_phi;

    // --- Device output buffers ---
    long Nsteps_total = (long)N * max_steps;
    StructureStep *d_bwd, *d_fwd;
    int *d_nsteps_bwd, *d_nsteps_fwd;
    CUDA_CHECK(cudaMalloc(&d_bwd,       Nsteps_total * sizeof(StructureStep)));
    CUDA_CHECK(cudaMalloc(&d_fwd,       Nsteps_total * sizeof(StructureStep)));
    CUDA_CHECK(cudaMalloc(&d_nsteps_bwd, N * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_nsteps_fwd, N * sizeof(int)));

    int threads = 256;
    int blocks  = (N + threads - 1) / threads;

    printf("GPU: launching structure kernels (%d points, %d max output steps) ...\n",
           N, max_steps);
    fflush(stdout);

    // Backward pass
    if(params.MapDirection <= 0) {
        k_structure<<<blocks, threads>>>(d_initR, d_initZ, d_initPhi,
                                          d_bwd, d_nsteps_bwd,
                                          N, max_steps,
                                          params.dpinit, params.nstep, params.itt,
                                          -1);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    } else {
        cudaMemset(d_nsteps_bwd, 0, N*sizeof(int));
    }

    // Forward pass
    if(params.MapDirection >= 0) {
        k_structure<<<blocks, threads>>>(d_initR, d_initZ, d_initPhi,
                                          d_fwd, d_nsteps_fwd,
                                          N, max_steps,
                                          params.dpinit, params.nstep, params.itt,
                                          +1);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    } else {
        cudaMemset(d_nsteps_fwd, 0, N*sizeof(int));
    }

    printf("GPU: structure kernels complete, retrieving results ...\n");
    fflush(stdout);

    // --- Copy back ---
    CUDA_CHECK(cudaMemcpy(results_bwd, d_bwd,       Nsteps_total*sizeof(StructureStep), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(results_fwd, d_fwd,       Nsteps_total*sizeof(StructureStep), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(steps_bwd,  d_nsteps_bwd, N*sizeof(int),                      cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(steps_fwd,  d_nsteps_fwd, N*sizeof(int),                      cudaMemcpyDeviceToHost));

    // --- Free device memory ---
    cudaFree(d_BR); cudaFree(d_Bphi); cudaFree(d_BZ); cudaFree(d_psi);
    if(d_wall_R) cudaFree(d_wall_R);
    if(d_wall_Z) cudaFree(d_wall_Z);
    cudaFree(d_initR); cudaFree(d_initZ); cudaFree(d_initPhi);
    cudaFree(d_bwd);   cudaFree(d_fwd);
    cudaFree(d_nsteps_bwd); cudaFree(d_nsteps_fwd);

    return 0;
}
