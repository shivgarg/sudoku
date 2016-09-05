Sudoku-Solver
============

- Sudoku solver written in C. The main approach is applying depth first search to look for solution while applying various heuristics which reduce the arity of each step.The solver uses openmp threads to explore several solution possibilites to have a speedup over traditional single threaded application.  Some of the heuristics applied :-
1.	Elimination: If a cell can only have a single number , then the cell is assigned that number.
2.	Lone Range (both row, column and sub-matrix): If in a given row ,column or sub-matrix a number can come in only one cell, then that cell is assigned the corresponding number. 
3.	Lone rangers in pairs, triplets: Extending the above logic to pair,triplets of numbers which helps to reduce arity for dfs.

### Compiling 
Run the command **make** . It produces an executable sudoku_solver. 

### Running
The executable takes file name as input which contains the borad status. '0' indicates an empty cell in sudoku matrix rest indicate the position is occupied.

./sudoku-solver num_threads filename


