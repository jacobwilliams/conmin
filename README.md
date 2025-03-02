This is a modern Fortran version of the CONMIN optimizer.

### Introduction

In many mathematical problems, it is necessary to determine the minimum or maximum of a function of several variables, limited by various linear and nonlinear inequality constraints. It is seldom possible, in practical applications, to solve these problems directly, and iterative methods are used to obtain the numerical solution. Machine-calculation of this solution is, of course, desirable and the CONMIN program has been developed to solve a wide variety of such problems.

CONMIN is a FORTRAN program, in subroutine form, for the minimization of a multi-variable function subject to a set of inequality constraints. The general minimization problem is: Find values for the set of variables, X(I), to

```
  Minimize    OBJ

  Subject to:

       G(J) <= 0                  J = 1,NCON

       VLB(I) <= X(I) <= VUB(I)   I = 1,NDV      NSIDE > 0
```

where OBJ is a general function (objective function) of the variables, X(I) referred to hereafter as decision variables. OBJ need not be a simple analytic function, and may be any function which can be numerically evaluated.

G(J) is the value of the Jth inequality constraint, which is also a function of the X(I). NCON is the number of constraints, G(J). NCON may be zero. VLB(I) and VUB(I) are lower and upper bounds respectively on variable X(I), and are referred to as side constraints. NSIDE = 0 indicates that no lower or upper bounds are prescribed. If NCON = 0 and NSIDE = 0, the objective function is said to be unconstrained. NDV is the total number of decision variables, X(I).

Constraint G(J) is defined as active if CT<=G(J)<=ABS(CT) and violated if G(J)>ABS(CT), where constraint thickness, CT, is a small specified negative number. The numerical significance of CT may be understood by referring to Fig. 1, which shows a single constraint in a two variable design space. Constraint G(J) is mathematically equal to zero along a single curve in design space. However, on a digital computer, the exact value of G(J) = 0 can seldom be obtained. Therefore, the "curve" becomes a thick band with constraint thickness of 2*ABS(CT) over which G(J) is assumed to be zero. Because all G(J) must be negative, CT is taken as a negative number for consistency so that any G(J)>CT is defined as active (or violated) if G(J)>ABS(CT). While it may seem logical to choose a very small value (in magnitude) for CT (say -1.0E-6), the nature of the optimization algorithm used by CONMIN is such that more numerical stability can be achieved by taking CT = -0.1 or even -0.2. CT is used for numerical stability only, and when the optimization is complete, one or more constraints, G(J), will usually be very near zero, as seen in the examples in SECTION VIII.

It is desirable that G(J) be normalized so that
```
  -1 <= G(J) <= 1             J = 1, NCON
```

In this way the constraint thickness, CT, has the same numerical significance for all G(J). It is not necessary that all G(J) be precisely in this form, and such normalization may not be possible. However, it is important that all G(J) at least be of the same order of magnitude. For example, assume that some G(J) = X(1)**2-X(2). If X(1) and X(2) are expected to be of order 100 for the particular problem under consideration, G(J) may be scaled by dividing by 10,000 to provide a value for G(J) of order one.

The basic analytic technique used by CONMIN is to minimize OBJ until one or more constraints, G(J), become active. The minimization process then continues by following the constraint boundaries in a direction such that the value of OBJ continues to decrease. When a point is reached such that no further decrease in OBJ can be obtained, the process is terminated. The value of the constraint thickness parameter, CT, and the normalization of the constraints, G(J), have considerable effect on the numerical stability and rate of convergence of the optimization process.

An example of a constrained nonlinear problem is the minimization of the four variable Rosen-Suzuki function (ref. 1):
```
     MINIMIZE OBJ = X(1)**2 - 5*X(1) + X(2)**2 - 5*X(2) +
                    2*X(3)**2 - 21*X(3) + X(4)**2 + 7*X(4) + 50

     Subject to:

         G(1) = X(1)**2 + X(1) + X(2)**2 - X(2) +
                X(3)**2 + X(3) + X(4)**2 - X(4) - 8  <= 0

         G(2) = X(1)**2 - X(1) + 2*X(2)**2 + X(3)**2 +
                2*X(4)**2 - X(4) - 10                <= 0

         G(3) = 2*X(1)**2 + 2*X(1) + X(2)**2 - X(2) +
                X(3)**2 - X(4) - 5                   <= 0
```

This problem has four decision variables and three constraints, (NDV = 4, NCON = 3). No lower or upper bounds VLB(I) or VUB(I) are prescribed so control parameter NSIDE is specified as NSIDE = 0 to indicate this. It is necessary to provide a set of initial values for X(I), and from this the constrained optimum is obtained by CONMIN and its associated routines.

The minimization algorithm is based on Zoutendijk's method of feasible directions (ref. 2). The algorithm has been modified to improve efficiency and numerical stability and to solve optimization problems in which one or more constraints, G(J), are initially violated (ref. 3). While the program is intended primarily for the efficient solution of constrained functions, unconstrained functions may also be minimized (NCON = 0 and NSIDE = 0), and the conjugate direction method of Fletcher and Reeves (ref. 4) is used for this purpose. If a function is to be maximized, this may be achieved by minimizing the negative of the function.

For constrained minimization problems, the initial design need not be feasible (one or more G(J) may be greater than ABS(CT)), and a feasible solution (if one exists) is obtained with a minimal increase in the value of the objective function.

The user must supply a main program to call subroutine CONMIN along with an external subroutine to evaluate the objective function, constraint functions and the analytic gradient of the objective and currently active or violated constraint functions. At any given time in the minimization process, gradient information is required only for constraints which are active or violated (G(J)>=CT). Gradients are calculated by finite difference if this information is not directly obtainable, and a subroutine is included with CONMIN for this purpose.


