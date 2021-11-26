# Finite Difference Method for the Multi-Asset Black–Scholes Equations
The code runs through the driver code in src/driver_code.m. The code has been implemented for 3 dimensions. The driver code asks for the dimension as input and returns the plots and error analysis. 

## Folder Structure:

```bash
│
├── outputs                             # contains the outputs for the executed code
│   ├── plots
│   │   └── ...
│   └── outputs.txt                 
│
├── papers                              # reference papers
│   └── ...
│
├── src
│   ├── cash_or_nothing_dimension1.m    # numerical solution for the cash of nothing option involving 1 asset
│   ├── cash_or_nothing_dimension2.m    # numerical solution for the cash of nothing option involving 2 assets
│   ├── cash_or_nothing_dimension3.m    # numerical solution for the cash of nothing option involving 3 assets
│   ├── closed_form_solution_1d.m       # closed form solution for the cash of nothing option involving 1 asset
│   ├── closed_form_solution_2d.m       # closed form solution for the cash of nothing option involving 2 assets
│   ├── closed_form_solution_3d.m       # closed form solution for the cash of nothing option involving 3 assets
│   ├── compute_error.m                 # calculates error between the numerical and the closed form solution
│   ├── driver_code.m                   # driver code
│   ├── plot_graphs.m                   # helper function to plot the figures
│   └── tridiagonal_matrix.m            # helper function to construct the tridiagonal matrix
│
└── term paper                          # term paper pdf and the presentation
    └── ... 
```