# Analysis-of-the-dynamics-of-a-transmission-tower-
Analysis of the dynamics of a transmission tower using MATLAB. 


by modelling the transmission tower in Matlab using 3D beam finite elements.
With each model, we are going to compute the first eight Eigen frequencies and draw the three first mode shapes.
Furthermore, we are going to model the damping in the
system using the proportional damping assumption to get a modal damping ratio of 1% for the
first two modes, then also listing the other damping ratio up to the eight one.


In the second part, we are going to analyse what happens when the power lines are
oscillated rapidly. To do this we will apply a vertical harmonic excitation characterised by
amplitude of 0.4KN and a period of 0.4s to the lumped masses. Then, we will compute
and compare the stationary response on a node: using modal
displacement method, acceleration method as well as the Newmark method.


Finally - modal reduction: we are going to reduce the model to compare the result with the initial
model. The purpose of this method is to lower the computation cost (time and storage) in Matlab
by reducing the size of the model.
To do this, we will decrease the number of degrees of freedom, by having a part of the structure
condensed and only the ten nodes retained. Using these nodes, we will build a
super-element with the methods: Guyan-Irons and Craig and Bampton’s such that the error on
the eight frequencies previously identified (with the initial model) is less than 0.2%.
