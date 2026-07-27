// GPU-compatible plain-C data structures for field-line tracing.
// No blitz++, no C++ stdlib — safe to include from both CUDA and g++ translation units.

#ifndef GPU_FIELDS_HXX
#define GPU_FIELDS_HXX

// ---------------------------------------------------------------------------
// 3-D B-field grid in cylindrical coordinates (R, phi, Z).
// Memory layout: index(iR, iphi, iZ) = iR*(Nphi*NZ) + iphi*NZ + iZ  (row-major)
// The psi grid is 2-D (axisymmetric):  index(iR, iZ) = iR*NZ + iZ
// All pointer fields are plain heap pointers — managed entirely on the CPU side;
// the CUDA wrappers copy them into device memory before launching any kernel.
// ---------------------------------------------------------------------------
struct FieldGrid3D
{
    // B-field components [NR * Nphi * NZ], Tesla
    double *BR;
    double *Bphi;
    double *BZ;

    // Normalised poloidal flux psi on 2-D (R,Z) grid [NR * NZ]
    // psi = (psi_raw - psi_axis) / (psi_sep - psi_axis)
    double *psi_norm;

    // Grid dimensions
    int NR, Nphi, NZ;

    // Grid extents
    double Rmin,   Rmax;
    double phimin, phimax;   // degrees, right-handed
    double Zmin,   Zmax;

    // Grid spacings (derived from extents and dims)
    double dR, dphi, dZ;

    // Wall polygon for boundary check [Nwall] each
    // wall_R[k] / wall_Z[k] are the k-th vertex, alternating order preserved from EQD.wall
    double *wall_R;
    double *wall_Z;
    int    Nwall;

    // Simple bounding box (from machine header bndy[])
    double bndy[4];   // [Rmin, Rmax, Zmin, Zmax]
    int    simpleBndy;

    // Magnetic axis (needed to compute theta on host)
    double RmAxis, ZmAxis;
};

// ---------------------------------------------------------------------------
// Initial condition for a single field line
// ---------------------------------------------------------------------------
struct FieldlineInit
{
    double R;    // metres
    double Z;    // metres
    double phi;  // degrees, right-handed
};

// ---------------------------------------------------------------------------
// Per-point scalar results returned by the laminar GPU kernel
// ---------------------------------------------------------------------------
struct LaminarResult
{
    double ntor;    // number of toroidal turns (both directions summed)
    double Lc;      // total connection length [m]
    double psimin;  // minimum normalised psi reached
    double psiav;   // mean normalised psi along trajectory
};

// ---------------------------------------------------------------------------
// Per-step results returned by the structure GPU kernel.
// Allocated as a flat [N * max_steps] array; kernel writes valid=0 after the
// field line has terminated so the host can detect the endpoint.
// ---------------------------------------------------------------------------
struct StructureStep
{
    double R;
    double Z;
    double phi;      // degrees, right-handed
    double psi;      // normalised
    double Lc;       // arc length from start [m]
    double dpsidLc;  // dpsi/dLc [1/m], signed with direction
    int    valid;    // 1 = step was computed, 0 = field line already terminated
};

// ---------------------------------------------------------------------------
// Runtime parameters passed to all GPU kernels
// ---------------------------------------------------------------------------
struct GPUTraceParams
{
    double dpinit;       // integration step size [degrees]
    int    nstep;        // structure only: output every nstep integration steps
    int    MapDirection; // +1 forward, -1 backward, 0 both
    int    itt;          // max toroidal iterations (laminar) or max output steps (structure)
    double phistart;     // starting toroidal angle [degrees]
    double v_par;        // prescribed parallel velocity [m/s] (drift mode; sets angle scale)
    double v_radial;     // prescribed anomalous radial velocity [m/s]; 0 -> off (pure field line)
    double v_tor;        // prescribed toroidal rotation velocity [m/s]; 0 -> off

    // Relativistic guiding-center drift (MAFOT particle mode).  sigma==0 -> off.  The constants
    // GAMMA/eps0/Ix are computed once on the host (PARTICLE) and replicate the CPU dgls drift.
    int    sigma;        // +1 co-passing, -1 counter-passing, 0 field lines only
    int    Zq;           // charge number (electrons -1, ions >=1)
    double GAMMA;        // relativistic gamma factor
    double eps0;         // normalized rest energy
    double Ix;           // normalized canonical angular momentum
    double R0;           // EFIT major radius R0 [m]
};

#endif // GPU_FIELDS_HXX
