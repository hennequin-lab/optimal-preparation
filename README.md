# optimal-preparation
Optimal anticipatory control as a theory of motor preparation: a thalamo-cortical circuit model

## Dependencies
1. [OCaml](https://ocaml.org) 4.10.0
2. [Owl](https://github.com/owlbarn/owl)
3. [Owl Opt](https://github.com/owlbarn/owl_opt)
4. [Cmdargs](https://github.com/hennequin-lab/cmdargs)

## Construct target reaches and ISN network
```sh
mkdir results
dune exec construct/reaches.exe -- -d results # create target reaches
dune exec soc/construct.exe -- -d results # create ISN
```
## Find target initial states and readout matrix
```sh
dune exec construct/setup.exe -- -d results # learn initial states and readout
```

## Feedback vs feedforward preparatory control
```sh
dune exec construct/vanilla.exe -- -d results # vanilla LQR preparation
dune exec construct/naive.exe -- -d results # naive feedforward preparation
```


## Full thalamo-cortical loop + preparation strategy
```sh
dune exec construct/dynamic_setup.exe -- -d results # setup full thalamocortical loop
dune exec construct/dynamic.exe -- -d results
```

## Plotting the results with gnuplot


To plot the activity of neuron `2` under optimal feedback (LQR) control for all 8 reaches, use the following command inside a gnuplot terminal within the `results` folder:
 
```gnuplot
# activity of neuron 2
plot for [i=1:8] "x_vanilla_loop_r1_".i u 2 w l lc black
```

Similarly, to plot the corresponding hand trajectories and propsective motor cost in log-scale, use the following command:

```gnuplot
# hand trajectories (x-y coordinates)
plot for [i=1:8] "hand_vanilla_loop_r1_".i u 1:3 w l lc black 
```

```gnuplot
# prospective motor cost
set log y
set xrange [0:1000]
plot for [i=1:8] "cost_vanilla_loop_r1_".i u 1 w l lc black
```

To plot similar quantities for feedforward control or full thalamo-cortical feedback control, substitute `vanilla` for `naive` or `dynamic` in the plotting commands above.
