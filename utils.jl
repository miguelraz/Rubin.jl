##################################################################
##  Chapter 1: Loading the damn thing.
#
# File for parsing the RubiRules and dumping them into a pretty printed JSON
# in your Rubin.jl/src dir.
#
# Here are the commands to download the artifacts:
# RUN THEM!
# using Pkg
# Pkg.instantiate()
##################################################################

using Rubin
using Test
using Artifacts
using FileTrees
using JSON3
using StructTypes
using Logging
using SymPy
using PyCall

# TODO mark as constants
rubi = artifact"Rubi"
rubitests = artifact"RubiTests"
rubidir = FileTree(rubi)
abstract type AbstractRubiParser end
struct RubiRules <: AbstractRubiParser end
struct RubiTests <: AbstractRubiParser end

## Design goals:
#1. Read Rules into a huge JSON array (with FileTrees.jl)
#2. Add metadata, facts, number, filename
#3. Save to a single file

## Steps:
# 1. Read artifact dir into a FileTree
# 2. Come up with a function to parse the files
# INPUT examples:
# (* ::Subsection::Closed:: *)
# (* 1.1.1.1 (a+b x)^m *)
# Int[1/x_, x_Symbol] := Log[x]
# Int[x_^m_., x_Symbol] := x^(m + 1)/(m + 1) /; FreeQ[m, x] && NeQ[m, -1]
# (Int[...]) := (...) /; (...)
# ^ lhs          ^rhs     ^facts
#
# BUT we can also have
# (* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(2*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)
# (* Int[...] := (...) ;/ (...) *)
# ^comment syntax around lhs := rhs ;/ facts *) closing comment syntax
# Which means we have 2 grammars:
# Case 1:
# (^Int.?) := (.+) /; (.+)
# Case 2:
# (* (Int.?) := (.+)(|\/; (.+)) 
#
## SOLUTION:
# With the help of Regex101.com:
# (Int.+) := (.+)(?: \/; \*\))
# ^glob lhs,
#             ^glob rhs
#                  ^glob facts, if any exist, ignore if it doesn't
# Let's now write some tests for the regex.
regextest = [
"""(* ::Subsection::Closed:: *)""",
"""(* 1.1.1.1 (a+b x)^m *)""",
"""Int[1/x_, x_Symbol] := Log[x]""",
"""(* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(3*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)""",
"""(* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(2*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)""",
"""(* Int[1/(a_+b_.*x_^5),x_Symbol] := With[{r=Numerator[Rt[a/b,5]],  s=Denominator[Rt[a/b,5]]}, r/(5*a)*Int[1/(r+s*x),x] + 2*r/(5*a)*Int[(r-1/4*(1-Sqrt[5])*s*x)/(r^2-1/2*(1-Sqrt[5])*r*s*x+s^ 2*x^2),x] + 2*r/(5*a)*Int[(r-1/4*(1+Sqrt[5])*s*x)/(r^2-1/2*(1+Sqrt[5])*r*s*x+s^ 2*x^2),x]] /; FreeQ[{a,b},x] && PosQ[a/b] *)""",
"""(* Int[1/Sqrt[a_+b_.*x_^3],x_Symbol] := With[{q=Rt[b/a,3]}, -Sqrt[2]*(1+Sqrt[3])*(1+Sqrt[3]+q*x)^2*Sqrt[(1+q^3*x^3)/(1+Sqrt[3]+ q*x)^4]/(3^(1/4)*q*Sqrt[a+b*x^3])* EllipticF[ArcSin[(-1+Sqrt[3]-q*x)/(1+Sqrt[3]+q*x)],-7-4*Sqrt[3]]]  /; FreeQ[{a,b},x] && PosQ[a] *)""",
"""(* Int[Sqrt[c_+d_.*x_^2]/((a_+b_.*x_^2)*Sqrt[e_+f_.*x_^2]),x_Symbol] :=   Sqrt[c+d*x^2]*Sqrt[c*(e+f*x^2)/(e*(c+d*x^2))]/(a*Rt[d/c,2]*Sqrt[e+f* x^2])* EllipticPi[1-b*c/(a*d),ArcTan[Rt[d/c,2]*x],1-c*f/(d*e)] /;  FreeQ[{a,b,c,d,e,f},x] && PosQ[d/c] *)"""];

## Where the 🌠 magic happens 🌠
intrulesregexfacts = r"(Int.+) := (.+) \/; (.+)(?: \*\))";
intrulesregexsimple = r"^(Int.+) := (.+)";
commenttest = """(* Int[Sqrt[a_.+b_.*x_]*(A_.+B_.*x_)/(Sqrt[c_.+d_.*x_]*Sqrt[e_.+f_.*x_ ]*Sqrt[g_.+h_.*x_]),x_Symbol] :=  B*Sqrt[a+b*x]*Sqrt[e+f*x]*Sqrt[g+h*x]/(f*h*Sqrt[c+d*x]) - B*(b*g-a*h)/(2*f*h)*Int[Sqrt[e+f*x]/(Sqrt[a+b*x]*Sqrt[c+d*x]*Sqrt[g+ h*x]),x] + B*(d*e-c*f)*(d*g-c*h)/(2*d*f*h)*Int[Sqrt[a+b*x]/((c+d*x)^(3/2)*Sqrt[ e+f*x]*Sqrt[g+h*x]),x] /; FreeQ[{a,b,c,d,e,f,g,h,A,B},x] &&  EqQ[2*A*d*f-B*(d*e+c*f),0] *)""";

## Regex tests
# Gotta make sure we are slurping up the things we care about, without too much cruft.
vregex = match.(intrulesregexsimple, regextest);
rulesregex = match.(intrulesregexfacts, regextest);
@test isnothing(vregex[1])
@test isnothing(vregex[2])
@test vregex[3].captures[1] == "Int[1/x_, x_Symbol]"
@test vregex[3].captures[2] == "Log[x]"
@test rulesregex[end].captures[3] == " FreeQ[{a,b,c,d,e,f},x] && PosQ[d/c]"

"""
# Structs for parsing
- Useful in conjunction with StructTypes.jl and dumping everything into a JSON
- Fields:
	pathname: String
	filename: String
	rulenumber: Int
	commented: Bool
	lhs: String
	rhs: String
	givens: String (not atomize)
"""
Base.@kwdef struct IntRuleCapture
	pathname::String = ""
	filename::String = ""
	rulenumber::Int = ""
	comment::Bool = true
	rhs::String = ""
	lhs::String = ""
	givens::String	 = ""
end

## Auxiliary functions for later
lhs(cap::RegexMatch) = cap.captures[1]
rhs(cap::RegexMatch) = cap.captures[2]
givens(cap::RegexMatch) = cap.captures[3]
iscommented(cap::RegexMatch) = startswith(cap.captures[1], '(')

# Input: A file
# Output: An array of structs like IntRuleCapture
# TODO: Figure out rule numbers
function Base.parse(file, ::Type{RubiRules})
	regex = r"(Int.+) := (.+) \/; (.+)( \*\))?"
	# Capture the lines that fit the regex, skip the empty ones, put it all into 
	# a big array. Gotta love Julia!
	caps = [i for i in match.(regex, readlines(file)) if !isnothing(i)]
	vlhs = lhs.(caps)
	vrhs = rhs.(caps)
	vgivens = givens.(caps)
	vcomments = iscommented.(caps)
	rubi = artifact"Rubi"
	path = relpath(file, rubi) # Fix this global?
	filename = splitpath(file)[end]
	# Assemble and return 📦 ➡
	[IntRuleCapture(pathname = path, 
					filename = filename,
					rulenumber = 0, 
					lhs = vlhs[i], 
					rhs = vrhs[i],
					givens = vgivens[i],
					comment = vcomments[i]) for i in 1:length(caps)]
end

## Testing nothing's borked with finding the files and making a TreeFile of them
intfiles = rubidir["Rubi-4.16.1.0/Rubi/IntegrationRules"];
parsedstructs = FileTrees.load(intfiles) do file
	parse(string(path(file)), RubiRules)
end;
vstructs = reducevalues(vcat, parsedstructs);
@test 7032 == length(vstructs)

"""
Loads the RubiRules into a vector.
See `IntRuleCapture` for the fields.
We use FileTrees.jl, so starting Julia with multiple threads should speed this up,
but it feels like it could be faster ... ?
"""
function load(::Type{RubiRules})
	rubi = artifact"Rubi"
	intfiles = joinpath(rubi, "Rubi-4.16.1.0", "Rubi", "IntegrationRules")
	files = FileTree(intfiles)
	parsedstructs = FileTrees.load(files) do file
		parse(string(path(file)), RubiRules)
	end
	vstructs = reducevalues(vcat, parsedstructs)
	7032 == length(vstructs) || error("Please alert @miguelraz, something has 💥")
	vstructs
end

## Writing to JSON! At last!
# Defining this straight from the JSON3 documentation
StructTypes.StructType(::Type{IntRuleCapture}) = StructTypes.Struct();
# For debugging
#json = JSON3.write(vstructs);
"""
Write `vstructs` into a `intrules.json` file in Rubin.jl/src, where `vstructs`
is defined as
```julia
vstructs = load(RubiRules)
```
"""
function write(vstructs::Vector{IntRuleCapture}, ::Type{RubiRules})
	# hat tip to Jacob Quinn and the #data gang
	# Find the Rubin root dir in an OS independent way
	targetpath = joinpath(pkgdir(Rubin), "src", "rubirules.json")
	# JSON3 way to write to a file with 💃 pretty printing 💃
	open(targetpath, "w") do f
		JSON3.pretty(f, JSON3.write(vstructs))
		println(f)
	end
end

##################################################################
# Chapter 2:
#        The friggin' tests
#
# GOALS: ⚽
#1. Read MathematicaSyntaxTestSuite into a huge JSON Array
#2. Add metadata, facts, number, filename, test number
#3. Save to a single file.
#
## Note: They claim 72944 tests total in the website, 
# but we must care about the ones that are commented
#
#                  💪 Currently at 72957 💪, 
#                          with 8        missing
#                            Wermer book missing
##################################################################

## Regex tests
# Tricky, gotta make sure we don't miss any of them!
inttests = 
["{Sqrt[2*x + 1], x, 1, (1/3)*(1 + 2*x)^(3/2)}",
"(* {Sqrt[2*x + 1], x, 1, (1/3)*(1 + 2*x)^(3/2)} *)"];
inttestregex = r"{(.+),(.+),(.+),(.+)}";
res = match.(inttestregex, inttests);
@test res[1].captures[1] == "Sqrt[2*x + 1]"
@test res[1].captures[2] == " x"
@test res[1].captures[3] == " 1"
@test res[1].captures[4] == " (1/3)*(1 + 2*x)^(3/2)"
@test res[2].captures[1] == "Sqrt[2*x + 1]"
@test res[2].captures[2] == " x"
@test res[2].captures[3] == " 1"
@test res[2].captures[4] == " (1/3)*(1 + 2*x)^(3/2)"

## Header regex
# We'll likely need some metadata to alleviate pain when updating the codebase.
# This should help mitigate some of those aches.
headerregex = r"\(\*(.+)\*\)";
@test match(headerregex, "(*Integrands of the form x^m PolyLog[n, a x^q]*)").captures[1] == "Integrands of the form x^m PolyLog[n, a x^q]"

"""
 Sample input: {Sqrt[2*x + 1], x,   1, (1/3)*(1 + 2*x)^(3/2)}
               ^query		    ^var ^steps      ^optimal answer
# NOTE: 
grep '(* {' | wc detects 284 commented cases at time of writing

# Fields:
	- filename
	- path
	- header
	- integrand
	- variable
	- steps
	- optimal
	- iscomment
"""
Base.@kwdef struct IntRuleTest
	filename::String =""
	path::String = ""
	integrand::String = ""
	variable::String = ""
	steps::String = "" # Yeah sure, this could be an Int but I'm 🏈 punting.
	optimal::String = ""
	iscomment::Bool = false
end

# header: string - to print testset and @async regions
# tests: Vec{IntRuleTest} - to write all the tests
Base.@kwdef struct IntRuleTestSection
	filename::String = ""
	path::String = ""
	header::String = ""
	tests::Vector{IntRuleTest} = IntRuleTest[]
end

"""
 1. Find first line that is not `skipme`, parse and push into the header
 2. Slurp the subsequent integrals, capture them into a vector of `IntRuleTest`s
 3. When you encounter a header that is different from the previous, put 
		all the previous IntRuleTests into a `IntTestSectionCapture` array
	Headers: There is much metadata like `(* n > 0 *)` and `(*Integrands ... )`
	We don't wanna throw that away if possible
# TODO: finer metadata with the headers in the file
So that testset nesting can happen
"""
function Base.parse(file, ::Type{RubiTests})
	# TODO add nesting to subsubsections
	# This regex globs up everything (*InsideTheirCommentSyntax*)
	headerregex = r"\(\*(.+)\*\)"
	# This regex globs up the 4 entries in {a,format,like,this}
	# TODO Fix those 5 missing cases in issue
	inttestregex = r"{(.+),(.+),(.+),(.+)?}"
	rubi = artifact"RubiTests"
	path = relpath(file)

	header = ""
	tests = IntRuleTest[]
	for line in readlines(file)
		# Skip the `(* ::Subsection...` and `(* ::Subsubsection...` lines
		isempty(line) && continue	
		startswith(line, "(* ::") && continue
		# If it's a test, push it
		if occursin(r"{", line) && occursin(r"}", line)
				# TODO Handle this better? We're losing information here
				# Perhaps adding a hierarchy to the section or subsection is good
			try
				m = match(inttestregex, line).captures
			catch
				@warn "inttestregex did not capture int test line"
				@warn line
				@warn file
				continue
			end
			m = match(inttestregex, line).captures

			# We need to know 7 things:
			# 1. integrand
			# 2. variable to integrate
			# 3. steps taken by rubi
			# 4. optimal answer
			# 5. if it's commented
			# 6. path
			# 7. filename
			integrand, variable, steps, optimal = (m...,)
			iscomment = startswith(line, "(* {")
			path = relpath(file, rubi) # Fix this global?
			filename = splitpath(file)[end]
			#@info "pushed to test"
			push!(tests, IntRuleTest(integrand = integrand, variable = variable,
									 path = path, filename =  filename,
									 steps = steps, optimal = optimal,
									 iscomment = iscomment))
		elseif startswith("(*", line) && endswith("*)", line)
		# 🆙 Update the header
			header = match(headerregex, line).captures[1]
		else 
			continue
		end
	end
	# Remember! Push the last section!
	tests
end

## Test the test parsing function
# Because that is one chonky boii
#testfile = "/home/mrg/.julia/artifacts/1148cba18dae2f8939af8bc542233a48cc42cf19/MathematicaSyntaxTestSuite-4.16.0/5 Inverse trig functions/5.5 Inverse secant/5.5.2 Inverse secant functions.m";
testfile = "parsetest.m"
vtests = parse(testfile, RubiTests)
@test length(vtests) == 11
@test vtests[3].iscomment
@test vtests[9].iscomment


# TODO swap out the loads
function findtreefolder(::Type{RubiRules})
	rubi = artifact"Rubi"
	intfiles = joinpath(rubi, "Rubi-4.16.1.0", "Rubi", "IntegrationRules")
	files = FileTree(intfiles)
end
function findtreefolder(::Type{RubiTests})
	rubitests = artifact"RubiTests"
	intfiles = joinpath(rubitests, "MathematicaSyntaxTestSuite-4.16.0")
	files = FileTree(intfiles)
end
"""
	Loads the RubiRules into a vector.
See `IntTestSections` for the fields.
"""
function load(::Type{RubiTests})
	rubitests = artifact"RubiTests"
	intfiles = joinpath(rubitests, "MathematicaSyntaxTestSuite-4.16.0")
	fs = FileTree(intfiles)
	# Wester problems has an irregular format, see issue # 3
	rm(fs, r"LICENSE|README.md|Wester Problems.m")
	@info fs
	parsedstructs = FileTrees.load(fs) do file
		parse(string(path(file)), RubiTests)
	end
	vtests = reducevalues(vcat, parsedstructs)
	72957 == length(vtests) || error("Please alert @miguelraz, something has 💥")
	vtests
end

vtests = load(RubiTests);
@test length(vtests) == 72957

# For JSON3 writing
StructTypes.StructType(::Type{IntRuleTest}) = StructTypes.Struct();

"""
Write `vtests` into a `inttests.json` file in Rubin.jl/src, where `vtests`
is defined as
```julia
vtests = load(RubiRules)
```
"""
function write(vtests::Vector{IntRuleTest}, ::Type{RubiTests})
	# Hat tip to Jacob Quinn 🎩 and the #data gang
	# Find the Rubin root dir in an OS independent way
	targetpath = joinpath(pkgdir(Rubin), "src", "rubitests.json")
	# JSON3 way to write to a file with 💃 pretty printing 💃
	open(targetpath, "w") do f
		JSON3.pretty(f, JSON3.write(vtests))
		println(f)
	end
end


##################################################################
# Chapter 3:
#        Woflram -> Julia
#
# GOALS: ⚽    
# 1. OUTPUT: tests/rubitests.jl with all the test files written in plain Julia
#	- Format shown below
# 2. OUTPUT: src/rubirules.jl with all the rules written in plain Julia
#   - Format: show below
#
# Approaches:
# 1. Make a small dict with Mathematica => Julia translations, 
#	THEN search and replace every string on both sides
# 2. use a while loop with the replacings, if it doesn't change go on.
# OOOPS
#  - All fun and games until E^x and E^(x+2) showed up - unclear how to deal 
#    with parens 💩
# 3. Use a library like RBNF.jl
# 4. DERP! See if someone else has solved it first!
# Sympy has a Mathematica parser - from there it should be much simple.
#  - 
# Options:
#
# 1. Symbolics.jl vs Metatheory.jl .... vs Rust egraphs?
# # Consider how to propagate assumptions
# https://julialang.zulipchat.com/#narrow/stream/236639-symbolic-programming/topic/Integration/near/232350734
##################################################################

### Notes:
# Once I have a symbolic sympy expression, I can just turn it into a Julia
# Julia string with `str(y)`:
sympy = pyimport("sympy")
str = py"'Sin[x]^2'";
mathparser = sympy.parsing.mathematica.mathematica
y = maths(str)
res = string(y) # 🚀
@test "sin(x)^2" == res

## Parser hijinks:
# Ugh - the sympy.mathematica parser does not accept Int[ ... :(
# No worries 😅, we can pass in a PyDict with our prefs and have it work!
function parsetosympy(s,d)
	s = pystring(PyObject(s))
    f = sympy.parsing.mathematica.mathematica
    f(s,d)
end

# note the *x to slurp all args
sympydictfixes = PyDict(Dict("Int[*x]" => "integrate(*x)",
							 "ProductLog[*x]" => "productlog(*x)",
							 "Gamma[*x]" => "gamma(*x)",
							 #"E^x" => "exp(x)",
							 #"E^(*x)" => "exp(*x)",
							 ))
@test "integrate(x,x)" == parsetosympy("Int[x,x]",sympydictfixes)
"integrate(x,x)"

## Parser hijinks:
# 1. Need to strip underscores and `x_Symbol` 
# 2.

### OOPS:
# It is clear that many of the "steps" field are being mangled horribly
steps = [length(v.steps) == 2 for v in vtests]
@test count(steps) == length(vtests)



