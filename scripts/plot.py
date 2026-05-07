from utilities.utilsPoint import plot_point_study
from utilities.utilsField import plot_field_study, E2_formula, E_inf_formula

# -------------------------------------
# 1. Plots related to Semilagrangian study
# -------------------------------------

tol_labels = {
    "normVel-tol-1e-3.csv": r"$tol=10^{-3}$",
    "normVel-tol-1e-2.csv": r"$tol=10^{-2}$",
    "normVel-tol-1e-1.csv": r"$tol=10^{-1}$",
    "normVel-tol-1.csv": r"$tol=10^{0}$",
    "normVel-tol-10.csv": r"$tol=10^{1}$",
}

cfl_center_labels = {
    "normVel-cfl-1e-2-center.csv": r"$CFL=10^{-2}$",
    "normVel-cfl-1e-1-center.csv": r"$CFL=10^{-1}$",
    "normVel-cfl-1-center.csv":    r"$CFL=1$",
    "normVel-cfl-10-center.csv":   r"$CFL=10$",
    "normVel-cfl-100-center.csv":  r"$CFL=100$",
    "normVel-cfl-1000-center.csv": r"$CFL=1000$",
}

cfl_topRight_labels = {
    "normVel-cfl-1e-2.csv": r"$CFL=10^{-2}$",
    "normVel-cfl-1e-1.csv": r"$CFL=10^{-1}$",
    "normVel-cfl-1.csv":    r"$CFL=1$",
    "normVel-cfl-10.csv":   r"$CFL=10$",
    "normVel-cfl-100.csv":  r"$CFL=100$",
    "normVel-cfl-1000.csv": r"$CFL=1000$",
}

plot_point_study(
    pattern      = "../results/outputs/SemiLagrange/uniform/tolerance/point/normVel-tol-*.csv",
    labels        = tol_labels,
    print_mapping = True,
    legend_title  = "Tolerance",
    ylabel       = r"$\|\mathbf{v}\|$ [m/s]",
    xlabel       = r"Time [s]",
    hlines       = [(1.0, r"$u_{ref}=1$")],
    xlim         = (0, 100),
    ylim         = (0.0, 1.3),
    xscale       = "linear",
    yscale       = "linear",
    save_path  = "../results/outputs/SemiLagrange/uniform/CFL/point/cfl-uniform-sl",
)

'''
plot_point_study(
    pattern       = "../results/outputs/SemiLagrange/uniform/CFL/point/center/normVel-cfl-*.csv",
    labels        = cfl_center_labels,
    print_mapping = True,
    legend_title  = "Courant number",
    ylabel        = r"$\|\mathbf{v}\|$ [m/s]",
    xlabel        = r"Time [s]",
    hlines        = [(1.0, r"$u_{ref}=1$")],
    xlim          = (0, 80),
    ylim          = (0.99, 1.015),
    xscale        = "linear",
    yscale        = "linear",
    save_path     = "../results/outputs/SemiLagrange/uniform/CFL/point/cfl-uniform-sl",
)
'''

'''
plot_point_study(
    pattern       = "../results/outputs/SemiLagrange/uniform/CFL/point/topRightEdge/normVel-cfl-*.csv",
    labels        = cfl_topRight_labels,
    print_mapping = True,
    legend_title  = "Courant number",
    ylabel        = r"$\|\mathbf{v}\|$ [m/s]",
    xlabel        = r"Time [s]",
    hlines        = [(1.0, r"$u_{ref}=1$")],
    xlim          = (0, 200),
    xscale        = "linear",
    yscale        = "linear",
    save_path     = "../results/outputs/SemiLagrange/uniform/CFL/point/cfl-uniform-sl",
)
'''


# 1.2 von-karman analysis

vk_ref_labels = {
    "normVel-sl-vk-ref-top.csv": r"upper point",
    "normVel-sl-vk-ref-bottom.csv": r"lower point",
}

plot_point_study(
    pattern       = "../results/outputs/SemiLagrange/von-karman/*.csv",
    print_mapping = True,
    labels=vk_ref_labels,
    ylabel        = r"$\|\mathbf{v}\|$ [m/s]",
    xlabel        = r"Time [s]",
    hlines        = [(1.0, r"$u_{inlet}=1$")],
    xlim          = (0, 500),
    xscale        = "linear",
    yscale        = "linear",
    save_path     = "../results/outputs/SemiLagrange/von-karman",
)