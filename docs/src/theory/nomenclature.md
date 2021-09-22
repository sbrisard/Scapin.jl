# Nomenclature

- ``d``: dimension of the physical space (typically ``d=2, 3``)
- ``\Omega``: ``d``-dimensional unit-cell
- ``L_1,\ldots, L_d``: dimensions of the unit-cell: ``\Omega=(0, L_1)\times(0, L_2)\times\cdots\times(0, L_d)``
- ``\lvert\Omega\rvert=L_1L_2\cdots L_d``: volume of the unit-cell
- ``\tuple{n}``: ``d``-dimensional tuple of integers ``\tuple{n}=(n_1, n_2, \ldots, n_d)``
- ``\tuple{N}=(N_1, N_2, \ldots, N_d)``: size of the simulation grid
- ``\lvert N\rvert=N_1N_2\cdots N_d``: total number of cells
- ``h_i=L_i/N_i``: size of the cells (``i=1, \ldots, d``)
- ``\cellindices=\{0, \ldots, N_1-1\}\times\{0, \ldots, N_2-1\}\times\{0, \ldots, N_3-1\}``: set of cell indices
- ``\Omega_{\tuple{p}}``: cells of the simulation grid (``\tuple{p}\in\cellindices``)
