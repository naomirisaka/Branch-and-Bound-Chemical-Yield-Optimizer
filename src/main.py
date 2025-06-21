import heapq
import re
import json
from typing import Dict, List, Set, Tuple, Optional
from datetime import datetime

class Reaction:
    def __init__(self, equation: str, time: float):
        self.equation = equation.strip()
        self.time = time  
        self.reactants, self.products = self.parse_equation(equation)

    def parse_equation(self, eq: str) -> Tuple[Dict[str, float], Dict[str, float]]:
        try:
            lhs, rhs = eq.split('->')
            reactants = {}
            products = {}
            
            # Parse left side (reactants)
            for r in lhs.strip().split('+'):
                coeff, compound = self.extract_compound_with_coefficient(r.strip())
                if compound:
                    reactants[compound] = coeff
            
            # Parse right side (products)
            for p in rhs.strip().split('+'):
                coeff, compound = self.extract_compound_with_coefficient(p.strip())
                if compound:
                    products[compound] = coeff
            
            return reactants, products
        except Exception as e:
            print(f"Error parsing equation '{eq}': {e}")
            return {}, {}

    def extract_compound_with_coefficient(self, term: str) -> Tuple[float, str]:
        match = re.match(r'(\d*\.?\d*)\s*([A-Za-z0-9()_]+)', term.strip())
        if match:
            coeff_str, compound = match.groups()
            coeff = float(coeff_str) if coeff_str else 1.0
            return coeff, compound
        return 1.0, term.strip()

    def calculate_yield(self, target_product: str, extent: float) -> float:
        if target_product in self.products:
            return self.products[target_product] * extent
        return 0.0

    def get_limiting_extent(self, available_moles: Dict[str, float]) -> float:
        if not self.reactants:
            return 0.0
        
        max_extent = float('inf')
        for reactant, coeff in self.reactants.items():
            if reactant in available_moles:
                available = available_moles[reactant]
                possible_extent = available / coeff
                max_extent = min(max_extent, possible_extent)
            else:
                return 0.0  # Missing reactant
        
        return max_extent if max_extent != float('inf') else 0.0

    def __str__(self):
        return f"{self.equation} | time={self.time} min"

class Node:
    def __init__(self, level: int, total_yield: float, total_time: float, 
                 available_moles: Dict[str, float], bound: float, taken: List[bool]):
        self.level = level
        self.total_yield = total_yield
        self.total_time = total_time  # Now in minutes
        self.available_moles = available_moles.copy()
        self.bound = bound
        self.taken = taken.copy()

    def __lt__(self, other):
        return self.bound > other.bound

# Periodic table data
ATOMIC_WEIGHTS = {
    'H': 1.008, 'He': 4.003, 'Li': 6.94, 'Be': 9.012, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95,
    'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.906, 'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.414,
    'In': 114.818, 'Sn': 118.710, 'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293
}
MOLECULAR_WEIGHTS = {}

def parse_chemical_formula(formula: str) -> Dict[str, int]:
    formula = expand_parentheses(formula)

    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    element_counts = {}
    for element, count_str in matches:
        count = int(count_str) if count_str else 1
        element_counts[element] = element_counts.get(element, 0) + count
    
    return element_counts

def expand_parentheses(formula: str) -> str:
    while '(' in formula:
        start = -1
        for i, char in enumerate(formula):
            if char == '(':
                start = i
            elif char == ')':
                end = i
                break
        else:
            break  
        
        inside = formula[start+1:end]
    
        multiplier_match = re.match(r'(\d+)', formula[end+1:])
        multiplier = int(multiplier_match.group(1)) if multiplier_match else 1
        multiplier_len = len(multiplier_match.group(1)) if multiplier_match else 0
        
        expanded = ''
        element_pattern = r'([A-Z][a-z]?)(\d*)'
        for element, count_str in re.findall(element_pattern, inside):
            count = int(count_str) if count_str else 1
            new_count = count * multiplier
            expanded += element + (str(new_count) if new_count > 1 else '')
        
        formula = formula[:start] + expanded + formula[end+1+multiplier_len:]
    
    return formula

def calculate_molecular_weight(formula: str) -> float:
    if formula in MOLECULAR_WEIGHTS:
        return MOLECULAR_WEIGHTS[formula]
    
    try:
        element_counts = parse_chemical_formula(formula)
        total_weight = 0.0
        
        for element, count in element_counts.items():
            if element in ATOMIC_WEIGHTS:
                total_weight += ATOMIC_WEIGHTS[element] * count
            else:
                print(f"Unknown element: {element}")
                while True:
                    try:
                        weight = float(input(f"Enter atomic weight for {element} (g/mol): "))
                        if weight > 0:
                            ATOMIC_WEIGHTS[element] = weight
                            total_weight += weight * count
                            break
                        print("Atomic weight must be positive!")
                    except ValueError:
                        print("Please enter a valid number!")
        
        MOLECULAR_WEIGHTS[formula] = total_weight
        print(f"Calculated molecular weight for {formula}: {total_weight:.3f} g/mol")
        return total_weight
        
    except Exception as e:
        print(f"Error calculating molecular weight for {formula}: {e}")
        while True:
            try:
                mw = float(input(f"Enter molecular weight for {formula} (g/mol): "))
                if mw > 0:
                    MOLECULAR_WEIGHTS[formula] = mw
                    return mw
                print("Molecular weight must be positive!")
            except ValueError:
                print("Please enter a valid number!")

def get_molecular_weight(compound: str) -> float:
    return calculate_molecular_weight(compound)

def convert_to_moles(amount: float, unit: str, compound: str) -> float:
    unit = unit.lower().strip()
    
    if unit in ['mol', 'moles']:
        return amount
    elif unit in ['g', 'gram', 'grams']:
        mw = get_molecular_weight(compound)
        return amount / mw
    elif unit in ['kg', 'kilogram', 'kilograms']:
        mw = get_molecular_weight(compound)
        return (amount * 1000) / mw
    else:
        print(f"Unknown unit: {unit}")
        print("Supported units: mol, g, kg")
        print("Assuming grams...")
        mw = get_molecular_weight(compound)
        return amount / mw

def load_from_file(filename: str) -> Tuple[List[Reaction], Dict[str, float], float, str]:
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        reactions = []
        compounds = {}
        max_time = 0
        target_product = ""
        
        section = None
        
        for line in lines:
            if line.startswith('#') or line == '':
                continue
            
            if line.upper() == "REACTIONS:":
                section = "reactions"
                continue
            elif line.upper() == "COMPOUNDS:":
                section = "compounds"
                continue
            elif line.upper() == "PARAMETERS:":
                section = "parameters"
                continue
            
            if section == "reactions":
                parts = line.split('|')
                if len(parts) == 2:
                    equation = parts[0].strip()
                    time = float(parts[1].strip())
                    reactions.append(Reaction(equation, time))
            
            elif section == "compounds":
                parts = line.split()
                if len(parts) == 3:
                    compound = parts[0]
                    amount = float(parts[1])
                    unit = parts[2]
                    moles = convert_to_moles(amount, unit, compound)
                    compounds[compound] = moles
            
            elif section == "parameters":
                if line.startswith("MAX_TIME:"):
                    max_time = float(line.split(':')[1].strip())
                elif line.startswith("TARGET:"):
                    target_product = line.split(':')[1].strip()
        
        return reactions, compounds, max_time, target_product
    
    except Exception as e:
        print(f"Error loading file: {e}")
        return [], {}, 0, ""

def save_solution_to_file(filename: str, max_yield: float, selected_reactions: List[str], 
                         stats: Dict, max_time: float, target_product: str, 
                         initial_data: dict = None):
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("CHEMICAL REACTION YIELD OPTIMIZATION SOLUTION\n")
            f.write("=" * 50 + "\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("PROBLEM PARAMETERS:\n")
            f.write(f"Target Product: {target_product}\n")
            f.write(f"Time Limit: {max_time:.2f} minutes\n\n")

            if initial_data:
                f.write("INITIAL CONDITIONS:\n")
                if 'reactions' in initial_data:
                    f.write(f"Available Reactions: {len(initial_data['reactions'])}\n")
                if 'compounds' in initial_data:
                    f.write("Initial Compounds:\n")
                    for compound, moles in initial_data['compounds'].items():
                        f.write(f"  {compound}: {moles:.3f} mol\n")
                f.write("\n")
            
            f.write("OPTIMIZATION RESULTS:\n")
            f.write(f"Maximum Yield: {max_yield:.3f} mol of {target_product}\n")
            f.write(f"Total Time Used: {stats['total_time']:.2f} minutes\n")
            f.write(f"Time Utilization: {(stats['total_time']/max_time)*100:.1f}%\n")
            f.write(f"Efficiency: {stats['efficiency']:.3f} mol/min\n")
            f.write(f"Nodes Explored: {stats['nodes_explored']}\n\n")
            
            if target_product in MOLECULAR_WEIGHTS:
                mass_yield = max_yield * MOLECULAR_WEIGHTS[target_product]
                f.write(f"Yield in grams: {mass_yield:.2f} g\n\n")
            
            f.write(f"OPTIMAL REACTION SEQUENCE ({len(selected_reactions)} reactions):\n")
            f.write("Format: Reaction equation (extent: X, time: Y min)\n")
            for i, reaction_info in enumerate(selected_reactions, 1):
                f.write(f"{i}. {reaction_info}\n")
            
            if not selected_reactions:
                f.write("No reactions selected (infeasible solution)\n")
            
            f.write("\n" + "=" * 50 + "\n")
            
        print(f"Solution saved to: {filename}")
        
    except Exception as e:
        print(f"Error saving solution to file: {e}")

def can_apply(reaction: Reaction, available_moles: Dict[str, float]) -> bool:
    return all(reactant in available_moles and available_moles[reactant] > 0 
               for reactant in reaction.reactants)

def calculate_bound(node: Node, reactions: List[Reaction], max_time: float, target_product: str) -> float:
    if node.total_time > max_time:
        return node.total_yield
    
    total_yield = node.total_yield
    time = node.total_time
    available_moles = node.available_moles.copy()
    
    for i in range(node.level + 1, len(reactions)):
        reaction = reactions[i]
        if can_apply(reaction, available_moles):
            
            extent = reaction.get_limiting_extent(available_moles)
            if extent > 0:
                actual_time = reaction.time * extent
                if time + actual_time <= max_time:
                    time += actual_time
                    potential_yield = reaction.calculate_yield(target_product, extent)
                    total_yield += potential_yield
                    
                    for reactant, coeff in reaction.reactants.items():
                        available_moles[reactant] -= coeff * extent
                    for product, coeff in reaction.products.items():
                        if product in available_moles:
                            available_moles[product] += coeff * extent
                        else:
                            available_moles[product] = coeff * extent
    
    return total_yield

def branch_and_bound(reactions: List[Reaction], max_time: float, 
                    initial_moles: Dict[str, float], target_product: str) -> Tuple[float, List[str], Dict]:
    if not reactions:
        return 0, [], {}
    
    def efficiency(r):
        max_extent = r.get_limiting_extent(initial_moles)
        potential_yield = r.calculate_yield(target_product, max_extent)
        actual_time = r.time * max_extent
        return potential_yield / (actual_time + 0.1)
    
    reactions.sort(key=efficiency, reverse=True)
    
    print(f"\nStarting optimization with {len(reactions)} reactions...")
    print(f"Target product: {target_product}")
    print(f"Time limit: {max_time} minutes")
    print(f"Sorted reactions by efficiency:")
    for i, r in enumerate(reactions[:5]):
        eff = efficiency(r)
        print(f"   {i+1}. {r.equation} (efficiency: {eff:.3f} mol/min)")
    if len(reactions) > 5:
        print(f"   ... and {len(reactions)-5} more")

    queue = []
    root = Node(-1, 0, 0, initial_moles, 0, [])
    root.bound = calculate_bound(root, reactions, max_time, target_product)
    heapq.heappush(queue, root)

    max_yield = 0
    best_combination = []
    nodes_explored = 0
    
    print(f"Initial upper bound: {root.bound:.3f} mol")

    while queue:
        node = heapq.heappop(queue)
        nodes_explored += 1
        
        if node.bound <= max_yield:
            continue
            
        if node.level == len(reactions) - 1:
            continue

        i = node.level + 1
        reaction = reactions[i]

        if can_apply(reaction, node.available_moles):
            extent = reaction.get_limiting_extent(node.available_moles)
            if extent > 0:
                actual_time = reaction.time * extent
                new_time = node.total_time + actual_time
                
                if new_time <= max_time:
                    yield_increase = reaction.calculate_yield(target_product, extent)
                    new_yield = node.total_yield + yield_increase
                    
                    new_moles = node.available_moles.copy()
                    for reactant, coeff in reaction.reactants.items():
                        new_moles[reactant] -= coeff * extent
                    for product, coeff in reaction.products.items():
                        if product in new_moles:
                            new_moles[product] += coeff * extent
                        else:
                            new_moles[product] = coeff * extent
                    
                    taken_node = Node(i, new_yield, new_time, new_moles, 0, node.taken + [True])
                    taken_node.bound = calculate_bound(taken_node, reactions, max_time, target_product)
                    
                    if new_yield > max_yield:
                        max_yield = new_yield
                        best_combination = taken_node.taken.copy()
                        print(f"New best yield: {max_yield:.3f} mol (Time: {new_time:.1f} min, Extent: {extent:.2f}, Actual time: {actual_time:.1f} min)")
                    
                    if taken_node.bound > max_yield:
                        heapq.heappush(queue, taken_node)

        not_taken_node = Node(i, node.total_yield, node.total_time, node.available_moles, 0, node.taken + [False])
        not_taken_node.bound = calculate_bound(not_taken_node, reactions, max_time, target_product)
        
        if not_taken_node.bound > max_yield:
            heapq.heappush(queue, not_taken_node)

    selected_reactions = []
    total_time = 0
    
    current_moles = initial_moles.copy()
    
    for i, taken in enumerate(best_combination):
        if taken and i < len(reactions):
            reaction = reactions[i]
            extent = reaction.get_limiting_extent(current_moles)
            actual_time = reaction.time * extent
            
            selected_reactions.append(f"{reaction.equation} (extent: {extent:.2f}, time: {actual_time:.1f} min)")
            total_time += actual_time
            
            for reactant, coeff in reaction.reactants.items():
                current_moles[reactant] -= coeff * extent
            for product, coeff in reaction.products.items():
                if product in current_moles:
                    current_moles[product] += coeff * extent
                else:
                    current_moles[product] = coeff * extent

    stats = {
        'nodes_explored': nodes_explored,
        'total_time': total_time,
        'efficiency': max_yield / max(total_time, 0.1)
    }

    return max_yield, selected_reactions, stats

def input_reactions() -> List[Reaction]:
    print("\n" + "="*60)
    print("ENTER CHEMICAL REACTIONS")
    print("="*60)
    print("Enter your chemical reactions. Type 'done' when finished.")
    print("Format: reactants -> products (with coefficients)")
    print("Example: (CH2)3O5 + 2 H2O -> 7 C6H12O6")
    print("Example: 2 H2 + O2 -> 2 H2O")
    
    reactions = []
    reaction_count = 0
    
    while True:
        reaction_count += 1
        print(f"\nReaction {reaction_count}:")
        
        equation = input("Enter reaction equation (or 'done' to finish): ").strip()
        if equation.lower() == 'done':
            break
        
        try:
            while True:
                try:
                    time = float(input("Reaction time in minutes (default 120.0): ") or "120.0")
                    if time > 0:
                        break
                    print("Reaction time must be positive")
                except ValueError:
                    print("Please enter a valid number")
            
            reaction = Reaction(equation, time)
            reactions.append(reaction)
            
            print(f"Added reaction: {reaction.equation}")
            print(f"   Reactants: {dict(reaction.reactants)}")
            print(f"   Products: {dict(reaction.products)}")
            print(f"   Time: {reaction.time} minutes")
            
        except Exception as e:
            print(f"Error adding reaction: {e}")
            print("Please try again with correct format")
            reaction_count -= 1
    
    if not reactions:
        print("No reactions entered. Adding example reaction...")
        example = Reaction("(CH2)3O5 + 2 H2O -> 7 C6H12O6", 120.0)
        reactions.append(example)
    
    print(f"\nTotal reactions added: {len(reactions)}")
    return reactions

def input_compounds(reactions: List[Reaction]) -> Dict[str, float]:
    print("\n" + "="*60)
    print("ENTER AVAILABLE COMPOUNDS")
    print("="*60)
    print("Enter the compounds you have available with amounts.")
    print("Type 'done' when finished.")
    
    compound_moles = {}
    
    all_compounds = set()
    for reaction in reactions:
        all_compounds.update(reaction.reactants.keys())
        all_compounds.update(reaction.products.keys())
    
    if all_compounds:
        print(f"\nCompounds mentioned in your reactions: {', '.join(sorted(all_compounds))}")
    
    while True:
        compound = input("\nEnter compound name (or 'done' to finish): ").strip()
        if compound.lower() == 'done':
            break
        
        if not compound:
            continue
        
        while True:
            try:
                amount = float(input(f"Amount of {compound}: "))
                if amount >= 0:
                    break
                print("Amount must be non-negative")
            except ValueError:
                print("Please enter a valid number")
        
        unit = input("Unit (g, kg, mol): ").strip()
        
        moles = convert_to_moles(amount, unit, compound)
        compound_moles[compound] = moles
        print(f"Added {moles:.3f} mol of {compound}")
    
    if not compound_moles:
        print("No compounds entered. Adding default compounds...")
        compound_moles = { # in moles
            "(CH2)3O5": 10.0, 
            "H2O": 20.0       
        }
    
    print(f"\nAvailable compounds (in moles):")
    for compound, moles in sorted(compound_moles.items()):
        print(f"   {compound}: {moles:.3f} mol")
    
    return compound_moles

def input_constraints(reactions: List[Reaction]) -> Tuple[float, str]:
    """Get time constraints and target product from user"""
    print("\n" + "="*60)
    print("ENTER OPTIMIZATION PARAMETERS")
    print("="*60)
    
    all_products = set()
    for reaction in reactions:
        all_products.update(reaction.products.keys())
    
    if all_products:
        print(f"Available products from your reactions: {', '.join(sorted(all_products))}")
    
    while True:
        target = input("Enter target product to maximize: ").strip()
        if target:
            break
        print("Please enter a target product!")
    
    while True:
        try:
            max_time = float(input("Maximum time limit (minutes): "))
            if max_time > 0:
                break
            print("Time must be positive!")
        except ValueError:
            print("Please enter a valid number!")
    
    print(f"\nTarget product: {target}")
    print(f"Time limit: {max_time} minutes")
    
    return max_time, target

def print_results(max_yield: float, selected_reactions: List[str], stats: Dict, 
                 max_time: float, target_product: str):
    print("\n" + "="*70)
    print("OPTIMIZATION RESULTS")
    print("="*70)
    
    if max_yield > 0:
        print(f"Maximum Yield of {target_product}: {max_yield:.3f} mol")
        print(f"Total Time Used: {stats['total_time']:.2f} / {max_time:.2f} minutes")
        print(f"Efficiency (mol/minute): {stats['efficiency']:.3f}")
        print(f"Nodes Explored: {stats['nodes_explored']}")
        
        print(f"\nOptimal Reactions Selected ({len(selected_reactions)}):")
        for i, eq in enumerate(selected_reactions, 1):
            print(f"   {i}. {eq}")
        
        time_util = (stats['total_time'] / max_time) * 100
        print(f"\nTime Utilization: {time_util:.1f}%")
        
        if target_product in MOLECULAR_WEIGHTS:
            mass_yield = max_yield * MOLECULAR_WEIGHTS[target_product]
            print(f"Yield in grams: {mass_yield:.2f} g")
        
    else:
        print("No feasible solution found within the given constraints.")
        print("Try increasing the time limit or check compound availability.")

def main():
    print("="*70)
    print("CHEMICAL REACTION YIELD OPTIMIZER")
    print("Branch and Bound Algorithm with Time Constraints")
    print("Time unit: MINUTES")
    print("="*70)
    
    # Ask for input method
    print("\nChoose input method:")
    print("1. Interactive input")
    print("2. Load from file")
    
    choice = input("Enter choice (1-2): ").strip()
    
    if choice == "2":
        filename = input("Enter input filename (e.g., input.txt): ").strip()
        reactions, compound_moles, max_time, target_product = load_from_file(filename)
        
        if not reactions:
            print("Failed to load from file. Switching to interactive mode.")
            reactions = input_reactions()
            compound_moles = input_compounds(reactions)
            max_time, target_product = input_constraints(reactions)
        else:
            print(f"Loaded {len(reactions)} reactions from file")
            print(f"Target: {target_product}, Time limit: {max_time} minutes")
            initial_data = {'reactions': reactions, 'compounds': compound_moles}
    else:
        reactions = input_reactions()
        compound_moles = input_compounds(reactions)
        max_time, target_product = input_constraints(reactions)
        initial_data = {'reactions': reactions, 'compounds': compound_moles}
    
    print(f"\nRunning Branch and Bound optimization...")
    max_yield, selected_reactions, stats = branch_and_bound(
        reactions, max_time, compound_moles, target_product
    )
    
    print_results(max_yield, selected_reactions, stats, max_time, target_product)
    
    save_choice = input("\nDo you want to save the solution to a file? (y/n): ").strip().lower()
    if save_choice in ['y', 'yes']:
        output_filename = input("Enter output filename (e.g., solution.txt): ").strip()
        if not output_filename:
            output_filename = f"solution_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        
        save_solution_to_file(output_filename, max_yield, selected_reactions, 
                             stats, max_time, target_product, initial_data)

if __name__ == "__main__":
    main()