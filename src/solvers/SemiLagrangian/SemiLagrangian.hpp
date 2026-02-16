#pragma once
#include "../../core/Fields.hpp"
#include "../../core/Parameters.hpp"
#include "../../core/OutputWriter.hpp"
#include <memory>

class SemiLagrangian {
public:
    // Constructor takes only parameters
    SemiLagrangian(const Parameters& params);
    
    // Destructor
    ~SemiLagrangian();
    
    // Main solver - runs the entire simulation
    void Run();
    
    // Single time step - advect then project
    void Step();
    
    // Access to fields for external use
    Fields2D& GetFields() { return *fields; }
    const Fields2D& GetFields() const { return *fields; }

private:
    // Simulation parameters
    const Parameters& params;
    int nx, ny;
    varType dx, dy, dt;
    varType density;
    
    // Fields (owned by this class)
    Fields2D* fields;
    
    // Output writers (optional based on params)
    std::unique_ptr<OutputWriter> uWriter;
    std::unique_ptr<OutputWriter> vWriter;
    std::unique_ptr<OutputWriter> pWriter;
    std::unique_ptr<OutputWriter> divWriter;
    
    // Initialize output writers based on parameters
    void InitializeOutputWriters();
    
    // Write output at current timestep
    void WriteOutput(int step) const;
    
    // Advection methods
    void Advect() const;
    void traceParticleU(int i, int j, varType& x, varType& y) const;
    void traceParticleV(int i, int j, varType& x, varType& y) const;
    [[nodiscard]] varType interpolateU(varType x, varType y) const;
    [[nodiscard]] varType interpolateV(varType x, varType y) const;
    void getVelocity(varType x, varType y, varType& u, varType& v) const;
    
    // Projection methods
    void MakeIncompressible();
    void solveJacobi(int maxIters, double tol);
    void updateVelocities();
};
