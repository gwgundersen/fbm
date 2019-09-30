
# Commands to apply a Gaussian process model to the simple regression problem.

gp-spec rlog.gp 1 1 100 / 0.05:0.5 0.05:0.5
model-spec rlog.gp real 0.05:0.5

data-spec rlog.gp 1 1 / rdata@1:100 . rdata@101:200 .

gp-gen rlog.gp fix 0.5 0.1
mc-spec rlog.gp heatbath hybrid 20:4 0.5
gp-mc rlog.gp 100
