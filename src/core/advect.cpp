
class SemiLagrangian {
public:
    SemiLagrangian(Fields2D& fields)
        : fields(fields), 
          nx(fields.nx), 
          ny(fields.ny),
          dx(fields.dx), 
          dy(fields.dy), 
          dt(fields.dt) {}
    
    // Advect velocity field using semi-Lagrangian method
    void Advect();
    
private:
    Fields2D& fields;
    size_t nx, ny;
    varType dx, dy, dt;
    
    // Trace particle backward in time from staggered grid position
    void traceParticleU(size_t i, size_t j, varType& x, varType& y);
    void traceParticleV(size_t i, size_t j, varType& x, varType& y);
    
    // Interpolate velocity at arbitrary position
    varType interpolateU(varType x, varType y);
    varType interpolateV(varType x, varType y);
    
    // Get velocity at arbitrary position (for RK2 tracing)
    void getVelocity(varType x, varType y, varType& u, varType& v);
};
