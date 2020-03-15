"""Tooltips for oncv GUI."""

def oncv_tip(key):
    """Tip for variable key."""
    return _ONCVTIPS[key]


_ONCVTIPS = dict(
atsym="Atomic symbol",
z="Atomic number",
nc="Number of core states",
nv="Number of valence states",
iexc="""\
Exchange-correlation functional:
   1-Wigner
   2-Hedin-Lundquist
   3-Perdew-Wang-Ceperly-Alder
   4-Perdew-Burke-Enzerhof""",
psfile="Format of pseudopotential file, psp8 for ABINIT, upf for PWSCF",
n="Principal quantum number",
l="Angular momentum",
f="Occupancy (MUST be >0)",
lmax="Maximum angular momentum for which psp is calculated (<=3)",
rc="Core radius for this l",
ep="""\
Energy at which psp is generated (eigenvalue inserted for occupied
state in reference configuration, positive energy must be specified
for barrier-confined "scattering" state for unoccupied l <=lmax.
A small positive energy is usually  good (0.1-0.25 Ha).""",
ncon="""\
Number of constraints for pseudo wave function to match all-electron wave function at rc,
value + ncon-1 derivatives, must be between 3 and 5 ("M" in the paper, Eq.(6))""",
nbas="Number of basis functions.  Must be between ncon+2 and ncon+5 (`N` in the paper, Eqs.(4-5))",
qcut="Wave vector defining `residual energy` in the RRKJ method (`q_c` in the paper, Eq.(1)",
lloc="""\
Angular momentum whose semi-local psp is taken as local.
lloc=4 denotes a smooth polynomial continuation of the all-electron potential.
If lloc<=lmax, remaining data are ignored, but must be there (zeros are OK).""",
lpopt="""\
Type of polynomial continuation for lloc=4. values 1-5  permitted.
    1) match 2 derivatives, r^0,2,4
    2) match 2 derivatives, r^0,4,6
    3) match 3 derivatives, r^0,4,5,6
    4) match 3 derivatives, r^0,4,6,8
    5) match 3 derivatives, r^0,2,4,6""",
rc5="Info not available",
dvloc0="""\
Shift of r=0 potential from basic continuation (Ha) for lloc=4
depends on lpopt above as follows, with x=(r/rc)
    1) dvloc0*(1-x^2)^3
    2) dvloc0*(1-x^4)^3
    3-5) dvloc0*(1-x^4)^4""",
nproj="Number of projectors, 1 or 2. Automatically set to 0 for l=lloc",
debl="""\
Energy added to basic psp  energy ep for 2nd projector, automatically reset to match 2nd projector with 2nd bound state
at this l when it is occupied (i.e., the psp is generated for a corresponding-l shallow core state)""",
icmod="""\
  0: no non-linear core correction charge.
  1: smooth monotonic polynomial model core charge fit at "matching" rc following reference 35.""",
fcfact="""\
Radius for above determined by  rho_core(r)=fcfact*rho_pseudo_valence(r).
Values 0.25-0.5 are usually good (look at plots)""",
epsh1="""\
Lower energy limit for `phase-shift-like` log-derivative plot,
should be below the nearest core level for this l to make sure
there are no ghosts, but -2.0 Ha usually is OK""",
epsh2="Upper energy limit, 2.0 usually good",
depsh="Energy mesh interval for plot, 0.02 usually good enough",
rcfact="rcfact scales the crossover radius to determine the range of the Teter function.",
rlmax="""\
Maximum radius for Abinit psp code 8 format output.
Must be greater than maximum rc (including lloc=4 rc), but also determines
range of diagnostic plots, so  ~2*rcmax is usually good""",
drl="""\
Mesh spacing of linear radial mesh for Abinit output.
0.02 is good for "softer" psps, 0.01 is probably better with 1st row, 3d's, or semi-core psps.""",
ncnf="""\
Number of test configurations (<=4).
The reference config is always run first as a consistency check.
Core always is the reference core configuration""",
)
#   nvcnf (repeated ncnf times) number of valence states in this test configuration
#   n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)
