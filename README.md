# Wave on Sierpinski

Simple python project simulating solution to wave equation on Sierpinski gasket according to my understanding of calculus on fractals.

---

## Mathematical Background

### Calculus on Fractals

On smooth domains (like the real line or a surface), calculus is defined using limits of differences over smaller and smaller distances. Fractals unlike manifolds lack coordinate charts and differential structure, so classical calculus breaks down.

In such a limited case we can still define quadratic form and use it to define a “graph Laplacian” as an analogue of the Laplace operator. This operator encodes how functions diffuse or vibrate over the fractal.

## Relation Between the Laplacian and Energy Quadratic Form

The Laplacian matrix $L$ encodes how function values vary between neighboring nodes on the fractal graph. It is closely related to the **energy quadratic form** defined on functions over the fractal.

For a function $u$ defined on the nodes, the energy quadratic form is:

$
\mathcal{E}(u) = \frac{1}{2} \sum_{i,j} A_{ij} (u_i - u_j)^2
$

where $A_{ij}$ is the adjacency matrix (1 if nodes $i$ and $j$ are connected, 0 otherwise).

This can be rewritten using the Laplacian matrix as:

$
\mathcal{E}(u) = \frac{1}{2} u^\top L u
$

Intuitively:

- $\mathcal{E}(u)$ measures the **“roughness”** or **“tension”** of $u$ over the fractal.
- Smaller energy means $u$ varies little between connected nodes (smoother function).
- Larger energy means bigger differences between neighbors (more oscillatory).

This quadratic form serves as the **potential energy** in the wave equation.

---

### Physical Interpretation

- The kinetic energy is given by $\frac{1}{2} \|\dot{u}\|^2$, representing motion.
- The potential energy $\frac{1}{2} u^\top L u$ represents stored energy due to the "elastic tension" of the function $u$ across the fractal.
- The total energy $E(t)$ combines both and is conserved in the ideal wave propagation.

This connection between the graph Laplacian and the energy quadratic form allows the Laplacian to act as the fractal analogue of the classical Laplace operator governing waves and diffusion.


---

## How We Build the Laplacian

To simulate waves on the Sierpinski gasket, we construct a graph approximation:

- Each node = one point on the fractal
- Edges connect neighboring points

Given an adjacency matrix $ A $ for the graph:

$
A_{ij} =
\begin{cases}
1 & \text{if nodes } i,j \text{ are connected}\\
0 & \text{otherwise}
\end{cases}
$

We define the **graph Laplacian**:

$
L = D - A
$

where $ D $ is the degree matrix:

$
D_{ii} = \text{number of neighbors of node } i
$

---

### Renormalization

Because the Sierpinski gasket shrinks distances at each level of its construction, the Laplacian also rescales. A key constant is the **renormalization factor**:

$
r = \frac{3}{5}
$

If we build the gasket to level $ \ell $, we multiply the Laplacian by:

$
\frac{1}{r^\ell}
$

This ensures the Laplacian approximates the continuum operator as the mesh becomes finer.

---

## The Wave Equation on a Fractal

We solve the wave equation:

$
\frac{\partial^2 u}{\partial t^2}
=
- L \, u
$

where:
- $ u $ = vector of values at each node
- $ L $ = Laplacian matrix of the fractal graph

This models how vibrations or waves propagate over the fractal.

---

## Spectral Solution

We diagonalize $ L $:

$
L \, \phi_k = \lambda_k \, \phi_k
$

Each eigenvector $ \phi_k $ describes a **mode** of vibration with frequency $ \sqrt{\lambda_k} $.

Given initial displacement $ f $ and velocity $ v $, we expand them:

$
f_k = \langle f, \phi_k \rangle
\quad
v_k = \langle v, \phi_k \rangle
$

Then the solution is:

$
u(t) = \sum_k
\left[
f_k \cos(\sqrt{\lambda_k}\, t)
+
\frac{v_k}{\sqrt{\lambda_k}} \sin(\sqrt{\lambda_k}\, t)
\right]
\phi_k
$

If $ \lambda_k = 0 $, the term becomes:

$
u_k(t) = f_k + v_k \cdot t
$

---

## Initial Conditions: Gaussian Bumps

Instead of exciting just one vertex, we use smooth **Gaussian bumps**:

$
f(x) =
A \cdot \exp
\left(
- \frac{\| x - x_{\text{source}} \|^2}{2\sigma^2}
\right)
$

Similarly for $v(x)$. This excites multiple eigenmodes, making waves visible.

---

## Energy

The total energy is:



$
E(t) =
\frac{1}{2} \left\| \dot{u}(t) \right\|^2
+
\frac{1}{2} \, u(t)^\top L \, u(t)
$


- Kinetic energy = $\frac{1}{2} | \dot{u} |^2$
- Potential energy = $\frac{1}{2} u^\top L u$

Energy should remain nearly constant over time (numerical errors may introduce tiny drifts).

---

## GIF Examples

### 2D Animation

![Wave propagation 2D](2d.gif)

---

### 3D Animation

![Wave propagation 3D](3d.gif)

---

## How to Run

1. Install dependencies:
    ```bash
    pip install numpy matplotlib scipy
    ```

2. Run the GUI:
    ```bash
    python sierpinski_wave_gui.py
    ```

