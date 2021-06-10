package engine

import "math"

//All variables here are global, because they are first letter capital.
var(
	R       [4]float64              // MTJ nanopillar resistance, Ohm, function of m(t)
	Rp      float64=float64(0.0)    // MTJ parallel resistance, Ohm, const
	Rap     float64=float64(0.0)    // SHO Resistance in AP state, Ohm, const
  Rout1   =NewScalarValue("Rout1","Ohm","SHNO 1 Resistance",func() float64 { return (Rap-Rp)*R[0]/2+(Rp+Rap)/2 }) //Use tableadd(Rout1) to print R, just as E_total
  Rout2   =NewScalarValue("Rout2","Ohm","SHNO 2 Resistance",func() float64 { return (Rap-Rp)*R[1]/2+(Rp+Rap)/2 }) //Note that R[n] is -mx (negative sign)
  Rout3   =NewScalarValue("Rout3","Ohm","SHNO 3 Resistance",func() float64 { return (Rap-Rp)*R[2]/2+(Rp+Rap)/2 })
  Rout4   =NewScalarValue("Rout4","Ohm","SHNO 4 Resistance",func() float64 { return (Rap-Rp)*R[3]/2+(Rp+Rap)/2 })
  Jread   float64=float64(0.0)    // MTJ read current density (A/m2)
  Gij   [4][4]float64           // Coupling conductance matrix
  Jcpl    =NewExcitation("Jcpl", "A/m2", "Coupling electrical current density") // Coupling current density added to mumax STT term
)

//Setting the variables using function calls in the script
func init() {
	DeclFunc("SetOutputResistance",SetOutputResistance,"Set MTJ P/AP resistance")
	DeclFunc("SetJread",SetJread,"Set MTJ read current density")
	DeclFunc("SetCoupling",SetCoupling,"Set coupling conductances between 4 SHNOs")
}

//Set MTJ parallel and anti-parallel resistances
func SetOutputResistance(rp_in float64, rap_in float64) {
  Rp = rp_in
  Rap = rap_in
}

//Set MTJ read current density
func SetJread(jread_in float64) {
  Jread = jread_in
}

// diagonal components of Gij matrix are 0, so oscillators do not couple with themselves
// Gij = Gij' (transpose)
// coupling is assumed reciprocal, s_ij = s_ji
func SetCoupling(s12 float64, s13 float64, s14 float64, s23 float64, s24 float64, s34 float64) {
  Gij[0][1] = s12
  Gij[0][2] = s13
  Gij[0][3] = s14
  Gij[1][0] = s12
  Gij[1][2] = s23
  Gij[1][3] = s24
  Gij[2][0] = s13
  Gij[2][1] = s23
  Gij[2][3] = s34
  Gij[3][0] = s14
  Gij[3][1] = s24
  Gij[3][2] = s34
}

// update MTJ resistance based on free layer magnetization (not yet scaled by P,AP resistances)
func UpdateOutputResistance() {
  for i:=0; i<4; i++ { // iterate over 4 SHNO
    if math.IsNaN(M.Region(i+1).Average().X()) { // if this region is not used in simulation ( < 4 SHNO)
        R[i] = float64(0.)
    } else {
        mx:=M.Region(i+1).Average().X()
        R[i]= -mx
    }
  }
  UpdateCouplingJ() // once MTJ resistance updated, use this to update coupling signal
}

// update coupling current density based on MTJ resistance as function of m
func UpdateCouplingJ() {
  Jcpl.SetRegionFn(1,func() [3]float64 { return [3]float64{0., 0., CalcCouplingJ(0)} })
  Jcpl.SetRegionFn(2,func() [3]float64 { return [3]float64{0., 0., CalcCouplingJ(1)} })
  Jcpl.SetRegionFn(3,func() [3]float64 { return [3]float64{0., 0., CalcCouplingJ(2)} })
  Jcpl.SetRegionFn(4,func() [3]float64 { return [3]float64{0., 0., CalcCouplingJ(3)} })
}

// calculate coupling current density for r region
func CalcCouplingJ(r int) float64 {
  gcpl := float64(0.0) // summed coupling signals
  for i:=0; i<4; i++ { // iterate over SHNOs coupled to r
    gcpl += ((Rap-Rp)*R[i]/2)*Gij[r][i] // do not include DC offset in coupling signal
  }
  return gcpl*Jread
}
