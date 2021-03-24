(* ::Package:: *)

(* ::Title:: *)
(*Integration Problems Involving Inverse Secants*)
(* ::Section::Closed:: *)
(*Integrands of the form x^m ArcSec[a x^n]*)

(* ::Subsection::Closed:: *)
(*Integrands of the form x^m ArcSec[a x^n]*)

(* ::Subsubsection::Closed:: *)
(*n>0*)

{ArcSec[a*x^5]/x, x, 7, (1/10)*I*ArcSec[a*x^5]^2 - (1/5)*ArcSec[a*x^5]*Log[1 + E^(2*I*ArcSec[a*x^5])] + (1/10)*I*PolyLog[2, -E^(2*I*ArcSec[a*x^5])]}
{x^3*ArcSec[Sqrt[x]], x, 4, (-(1/4))*Sqrt[-1 + x] - (1/4)*(-1 + x)^(3/2) - (3/20)*(-1 + x)^(5/2) - (1/28)*(-1 + x)^(7/2) + (1/4)*x^4*ArcSec[Sqrt[x]]}
(* {x^2*ArcSec[Sqrt[x]], x, 4, (-(1/3))*Sqrt[-1 + x] - (2/9)*(-1 + x)^(3/2) - (1/15)*(-1 + x)^(5/2) + (1/3)*x^3*ArcSec[Sqrt[x]]} *)
{ArcSec[Sqrt[x]]/x^4, x, 7, Sqrt[-1 + x]/(18*x^3) + (5*Sqrt[-1 + x])/(72*x^2) + (5*Sqrt[-1 + x])/(48*x) - ArcSec[Sqrt[x]]/(3*x^3) + (5/48)*ArcTan[Sqrt[-1 + x]]}

(* ::Subsubsection::Closed:: *)
(*n<0*)

{x^2*ArcSec[a/x], x, 5, (-(1/3))*a^3*Sqrt[1 - x^2/a^2] + (1/9)*a^3*(1 - x^2/a^2)^(3/2) + (1/3)*x^3*ArcCos[x/a]}
{x^1*ArcSec[a/x], x, 4, (-(1/4))*a*x*Sqrt[1 - x^2/a^2] + (1/2)*x^2*ArcCos[x/a] + (1/4)*a^2*ArcSin[x/a]}
{x^0*ArcSec[a/x], x, 3, (-a)*Sqrt[1 - x^2/a^2] + x*ArcCos[x/a]}
{ArcSec[a/x]/x^1, x, 6, (-(1/2))*I*ArcCos[x/a]^2 + ArcCos[x/a]*Log[1 + E^(2*I*ArcCos[x/a])] - (1/2)*I*PolyLog[2, -E^(2*I*ArcCos[x/a])]}
(* {ArcSec[a/x]/x^2, x, 5, -(ArcCos[x/a]/x) + ArcTanh[Sqrt[1 - x^2/a^2]]/a} *)
{ArcSec[a/x]/x^3, x, 3, Sqrt[1 - x^2/a^2]/(2*a*x) - ArcCos[x/a]/(2*x^2)}
{ArcSec[a/x]/x^4, x, 6, Sqrt[1 - x^2/a^2]/(6*a*x^2) - ArcCos[x/a]/(3*x^3) + ArcTanh[Sqrt[1 - x^2/a^2]]/(6*a^3)}

(* ::Subsection::Closed:: *)
