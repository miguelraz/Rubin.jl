### TODO - setup a proper build.jl script and not have this script run oon a local git clone of
# TODO extract rules form symja?

#https://github.com/RuleBasedIntegration/Rubi
#julia> add_artifact!("Artifacts.toml", "Rubi", "https://github.com/RuleBasedIntegration/Rubi/archive/4.16.1.0.tar.gz", force=true)

#https://github.com/RuleBasedIntegration/MathematicaSyntaxTestSuite
#
#julia> add_artifact!("Artifacts.toml", "RubiTests", "https://github.com/RuleBasedIntegration/MathematicaSyntaxTestSuite/archive/4.16.0.tar.gz", force=true)
#julia> Pkg.instantiate()
using Rubin
using Artifacts
using FileTrees
using JSON3
using Test

# How to test artifact is properly installed?
#@test run(`ls $(artifact"Rubi")`)
#Pkg.Artifacts.ensure_artifact_installed("Rubi", "Artifacts.toml")
#OR
#Pkg.Artifacts.artifact_exists(hash::SHA1)
#@test run(`ls $(artifact"RubiTests")`)
rubi = artifact"Rubi"
rubitests = artifact"RubiTests"
run(`ls $rubi`)

#0 Design goals:
#1. read Rules into a huge JSON array (with FileTrees.jl)
#2. Add metadata, facts, number, filename
#3. save to a single file

# Steps:
# 1. Read artifact dir into a FileTree
rubidir = FileTree(rubi)

# 2. Come up with a function to parse the files
# INPUT examples:
#(* ::Subsection::Closed:: *)
#(* 1.1.1.1 (a+b x)^m *)
#Int[1/x_, x_Symbol] := Log[x]
#Int[x_^m_., x_Symbol] := x^(m + 1)/(m + 1) /; FreeQ[m, x] && NeQ[m, -1]
# (Int[...]) := (...) /; (...)
# ^ lhs          ^rhs     ^facts
# BUT we can also have
# ...
# (* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(2*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)
# (* Int[...] := (...) ;/ (...) *)
# ^comment syntax around lhs := rhs ;/ facts *) closing comment syntax
# Which means we have 2 grammars:
# Case 1:
# (^Int.?) := (.+) /; (.+)
# Case 2:
# (* (Int.?) := (.+)(|\/; (.+)) 
#
#
# SOLUTION:
# With the help of Regex101.com:
# (Int.+) := (.+)(|\/; (.+))
# ^glob lhs,
#             ^glob rhs
#                  ^glob facts, if any exist
# REMEMBER: Add a marker if test is commented

regextest = [
"""(* ::Subsection::Closed:: *)""",
"""(* 1.1.1.1 (a+b x)^m *)""",
"""Int[1/x_, x_Symbol] := Log[x]""",
"""(* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(2*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)""",
"""(* Int[1/(a_+b_.*x_^5),x_Symbol] := With[{r=Numerator[Rt[a/b,5]],  s=Denominator[Rt[a/b,5]]}, r/(5*a)*Int[1/(r+s*x),x] + 2*r/(5*a)*Int[(r-1/4*(1-Sqrt[5])*s*x)/(r^2-1/2*(1-Sqrt[5])*r*s*x+s^ 2*x^2),x] + 2*r/(5*a)*Int[(r-1/4*(1+Sqrt[5])*s*x)/(r^2-1/2*(1+Sqrt[5])*r*s*x+s^ 2*x^2),x]] /; FreeQ[{a,b},x] && PosQ[a/b] *)""",
"""(* Int[1/Sqrt[a_+b_.*x_^3],x_Symbol] := With[{q=Rt[b/a,3]}, -Sqrt[2]*(1+Sqrt[3])*(1+Sqrt[3]+q*x)^2*Sqrt[(1+q^3*x^3)/(1+Sqrt[3]+ q*x)^4]/(3^(1/4)*q*Sqrt[a+b*x^3])* EllipticF[ArcSin[(-1+Sqrt[3]-q*x)/(1+Sqrt[3]+q*x)],-7-4*Sqrt[3]]]  /; FreeQ[{a,b},x] && PosQ[a] *)""",
"""(* Int[Sqrt[c_+d_.*x_^2]/((a_+b_.*x_^2)*Sqrt[e_+f_.*x_^2]),x_Symbol] :=   Sqrt[c+d*x^2]*Sqrt[c*(e+f*x^2)/(e*(c+d*x^2))]/(a*Rt[d/c,2]*Sqrt[e+f* x^2])* EllipticPi[1-b*c/(a*d),ArcTan[Rt[d/c,2]*x],1-c*f/(d*e)] /;  FreeQ[{a,b,c,d,e,f},x] && PosQ[d/c] *)"""]
intrulesregexfacts = r"(Int.+) := (.+)(?:\/; (.+) \*\))"
inturulesregexsimple = r"^(Int.+) := (.+)"

vregex = match.(intrulesregexsimple, regextest)
@test isnothing(vregex[1])
@test isnothing(vregex[2])
@test vregex[3].captures[1] == "Int[1/x_, x_Symbol]"
@test vregex[3].captures[2] == "Log[x]"
@test vregex[3].captures[3] == ""

function intrulesfileparser(file)
	
end

# 3. Lazy load the files
lazyjsons = FileTrees.load(rubidir; lazy = true) do file

end

#1. read MathematicaSyntaxTEstSuite into a huge JSON Array
#2. Add metadata, facts, number, filename, test number
#3. save to a single file.
#
#AFTER the parse and dump:
#Construct Mathematica to Julia rewriter?
#
# builtin reading/writing
# JSON3.read(json_string)
# JSON3.write(x)

# custom types
# JSON3.read(json_string, T; kw...)
# JSON3.write(x)
# builtin reading/writing
# JSON3.read(json_string)
# JSON3.write(x)
#
# # custom types
# JSON3.read(json_string, T; kw...)
# JSON3.write(x)
#
# More complicated
#
# # custom types: incrementally update a mutable struct
# x = T()
# JSON3.read!(json_string, x; kw...)
# JSON3.write(x)
#
# # read from file
# json_string = read("my.json", String)
# JSON3.read(json_string)
# JSON3.read(json_string, T; kw...)
#
# # write to file
# open("my.json", "w") do f
#     JSON3.write(f, x)
#         println(f)
#         end
#
#         # write a pretty file
#         open("my.json", "w") do f
#             JSON3.pretty(f, JSON3.write(x))
#                 println(f)
#                 end


# custom types: incrementally update a mutable struct
# x = T()
# JSON3.read!(json_string, x; kw...)
# JSON3.write(x)

# read from file
# json_string = read("my.json", String)
# JSON3.read(json_string)
# JSON3.read(json_string, T; kw...)

# write to file
# open("my.json", "w") do f
    # JSON3.write(f, x)
    # println(f)
# end

# write a pretty file
# open("my.json", "w") do f
    # JSON3.pretty(f, JSON3.write(x))
    # println(f)
# end
