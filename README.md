# Chemical Reaction Yields Optimization Using The Branch and Bound Algorithm
> Makalah IF2211 Strategi Algoritma

## Project Overview
<p align = justify>This project applies the Branch and Bound algorithm to optimize chemical reaction yields. Each node in the state-space tree represents a partial solution where decisions have been made for a subset of reactions, while each edge represents a decision on whether or not to execute a reaction. Time constraints, available compounds, and stoichiometric balance act as limiting factors, and a bounding function enables pruning of suboptimal paths. The program processes the input data, ranks reactions by effectiveness, and selects the most optimal sequence to maximize yield based on the defined constraints.</p>

## Project Structure
The project is structured as below:
```
├── .gitignore
├── README.md
├── src/
│   ├── io_handler.py     # Handles input/output
│   ├── main.py           # Main program
│   ├── optimizer.py      # Branch and Bound logic
│   └── reactions.py      # Data models and utilities
└── test/
```

## How to Use
1. Clone the repository: 
   ```sh
   git clone https://github.com/naomirisaka/Branch-and-Bound-Chemical-Yield-Optimizer.git

2. Open the cloned local repository in an IDE that supports Python.

3. Run the program:
   ```sh
   cd src
   python main.py
