#run this script to host the documentation locally
import Pkg ; Pkg.activate("docs")
include("docs/make.jl")
using LiveServer
serve(dir="docs/build")
