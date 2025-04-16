# Really important: This assumes balanced panels
using CSV
using DataFrames
using ProjectRoot
here = ProjectRoot.@projectroot
include("$here/code/qld/QLD.jl")
include("$here/code/qld/gmm_qld.jl")
include("$here/code/qld/attgt.jl")
include("$here/code/qld/within_transform.jl")
include("$here/code/qld/qld_imputation.jl")

df = DataFrame(CSV.File("data/Simulations/df_ex_dgp_3.csv"))
y = :y
id = :id
t = :t
g = :g
W = [:W1, :W2]
do_within_transform = false
# p = Int(2)
# p = Int(1)
p = Int(-1)
type = "dynamic"
return_y0 = false
return_naive_se = false


ret = qld_imputation(
  df;
  y=:y,
  id=:id,
  t=:t,
  g=:g,
  W=[:W1, :W2],
  do_within_transform=false,
  p=Int(-1),
  type="dynamic",
  return_y0=false,
  return_naive_se=false,
)
ret[:rel_year]
ret[:estimate]
sqrt.(diag(ret[:vcov]))
