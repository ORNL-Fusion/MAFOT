#include <cmath>
#include <exception>
#include <memory>
#include <string>

#include <mafot.hxx>

#if defined(_WIN32)
#define MAFOT_BRIDGE_EXPORT __declspec(dllexport)
#else
#define MAFOT_BRIDGE_EXPORT __attribute__((visibility("default")))
#endif

namespace {

thread_local std::string last_error;

struct MAFOTEngineHandle {
    EFIT EQD;
    IO PAR;
    COLLISION COL;
    std::unique_ptr<PARTICLE> FLT;

    MAFOTEngineHandle(
        const LA_STRING& run_directory,
        const LA_STRING& TprofileFile, const LA_STRING& NprofileFile,
        double f, double zbar, double rc, double mc,
        const char* boundary1, const char* boundary2,
        const bool collisions
    ) : PAR(EQD, run_directory + "/_iter.dat", 0) {
        double Raxis = 0, Zaxis = 0;
#ifdef m3dc1
        if (PAR.response_field == 0 || PAR.response_field == 2) {
            M3D.read_m3dc1sup(run_directory + "/");
            M3D.open_source(PAR.response, PAR.response_field, -1);
            Raxis = M3D.RmAxis;
            Zaxis = M3D.ZmAxis;
            M3D.unload();
        }
#endif

        EQD.ReadData(EQD.gFile, Raxis, Zaxis);

        if (boundary1 && *boundary1) EQD.setCustomBndy1(boundary1);
        if (boundary2 && *boundary2) EQD.setCustomBndy2(boundary2);

        if (collisions) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);

        prepare_common_perturbations(EQD, PAR, 0);
        prep_perturbation(EQD, PAR, 0, run_directory + "/");
        FLT.reset(new PARTICLE(EQD, PAR, COL));
    }

    ~MAFOTEngineHandle() noexcept {
#ifdef m3dc1
        if (PAR.response_field >= 0) M3D.unload();
#endif
    }
};

}  // namespace

extern "C" {

MAFOT_BRIDGE_EXPORT void* mafot_engine_create(
    const char* run_directory,
    const char* TprofileFile, const char* NprofileFile,
    double f, double zbar, double rc, double mc,
    const char* boundary1, const char* boundary2,
    const bool collisions
) noexcept {
    last_error.clear();
    if (!run_directory || !*run_directory) {
        last_error = "run_directory must not be null or empty";
        return nullptr;
    }

    try {
        return new MAFOTEngineHandle(
            LA_STRING(run_directory), LA_STRING(TprofileFile), LA_STRING(NprofileFile),
            f, zbar, rc, mc, boundary1, boundary2, collisions
        );
    } catch (const std::exception& exc) {
        last_error = exc.what();
    } catch (...) {
        last_error = "unknown error while creating MAFOTEngineHandle";
    }
    return nullptr;
}

MAFOT_BRIDGE_EXPORT int mafot_engine_map_once(
    void* handle, double r, double z, int direction, double* out
) noexcept {
    last_error.clear();
    if (!handle || !out) {
        last_error = "handle and output pointers must not be null";
    } else if (direction < -1 || direction > 1) {
        last_error = "direction must be -1, 0, or 1";
    } else if (!std::isfinite(r) || !std::isfinite(z)) {
        last_error = "R and Z must be finite";
    } else {
        auto* engine = static_cast<MAFOTEngineHandle*>(handle);
        if (direction == 0) {
            engine->FLT->set(1, 1, r, r, z, z, 1, 2);
            out[0] = engine->FLT->R;
            out[1] = engine->FLT->Z;
            return 0;
        }

        engine->FLT->R = r;
        engine->FLT->Z = z;
        engine->FLT->phi = engine->PAR.phistart;

        int chk = engine->FLT->mapit(1, direction, true);
        if (chk < 0) return chk;

        out[0] = engine->FLT->R;
        out[1] = engine->FLT->Z;
        engine->FLT->get_psi(engine->FLT->R, engine->FLT->Z, engine->FLT->psi);
        out[2] = engine->FLT->psi;
        out[3] = engine->FLT->get_theta();
        out[4] = engine->FLT->get_r();
        out[5] = engine->FLT->phi;
        return 0;
    }
    return -1;
}

MAFOT_BRIDGE_EXPORT void mafot_engine_destroy(void* opaque_handle) noexcept {
    delete static_cast<MAFOTEngineHandle*>(opaque_handle);
}

MAFOT_BRIDGE_EXPORT const char* mafot_engine_last_error() noexcept {
    return last_error.c_str();
}

}  // extern "C"
