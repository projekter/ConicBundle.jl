using Documenter
import ConicBundle

makedocs(sitename="ConicBundle.jl", modules=[ConicBundle], format=Documenter.HTML(prettyurls=false), expandfirst=["cref.md"])

deploydocs(repo="github.com/projekter/ConicBundle.jl.git")