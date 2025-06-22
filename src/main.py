from typing import Dict, List, Tuple
from reactions import Reaction, convert_to_moles, MOLAR_MASS
from optimizer import branch_and_bound

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
            
            if target_product in MOLAR_MASS:
                mass_yield = max_yield * MOLAR_MASS[target_product]
                f.write(f"Yield in grams: {mass_yield:.2f} g\n\n")
            
            f.write(f"OPTIMAL REACTION SEQUENCE ({len(selected_reactions)} reactions):\n")
            for i, reaction_info in enumerate(selected_reactions, 1):
                f.write(f"{i}. {reaction_info}\n")
            
            if not selected_reactions:
                f.write("No reactions selected (infeasible solution)\n")
            
        print(f"Solution saved to: {filename}")
        
    except Exception as e:
        print(f"Error saving solution to file: {e}")

def generate_input_example() -> str:
    try:
        with open('../test/example.txt', 'w', encoding='utf-8') as f:
            f.write("REACTIONS:\n")
            f.write("2 H2 + O2 -> 2 H2O | 5\n")
            f.write("C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O | 10\n")
            f.write("CuSO4 + 5 H2O -> CuSO4_5H2O | 3.5\n")
            f.write("\nCOMPOUNDS:\n")
            f.write("H2 10 mol\n")
            f.write("O2 5 mol\n")
            f.write("C6H12O6 1 mol\n")
            f.write("\nPARAMETERS:\n")
            f.write("MAX_TIME: 30\n")
            f.write("TARGET: H2O")
    except Exception as e:
        print(f"Error generating example input: {e}")

def input_reactions() -> List[Reaction]:
    print("\n" + "="*70 + "\n")
    print("ENTER CHEMICAL REACTIONS")
    print("Enter your chemical reactions. Type 'done' when finished.")
    print("Format: reactants -> products (with coefficients)")
    print("Example: (CH2)3O5 + 2 H2O -> 7 C6H12O6")
    
    reactions = []
    reaction_count = 0
    
    while True:
        reaction_count += 1
        print(f"\nReaction {reaction_count}:")
        
        while True:
            equation = input("Enter reaction equation (or 'done' to finish): ").strip()
            if equation.lower() == 'done':
                if reactions:  
                    return reactions
                else:
                    print("You must enter at least one reaction before finishing.")
                    continue
            
            if not equation:
                print("Please enter a valid reaction equation\n")
                continue
            
            time = None
            while time is None:
                time_input = input("Reaction time in minutes: ").strip()
                if not time_input:
                    print("Please enter a time value")
                    continue
                try:
                    time = float(time_input)
                    if time <= 0:
                        print("Reaction time must be positive")
                        time = None
                        continue
                except ValueError:
                    print("Please enter a valid number")
                    continue
            
            try:
                reaction = Reaction(equation, time)
                
                if not reaction.reactants and not reaction.products:
                    print("Error: Could not parse reaction equation. Please check the format.")
                    print("Make sure to use '->' to separate reactants and products\n")
                    continue
                
                if not reaction.reactants:
                    print("Error: No reactants found. Please check the equation format.")
                    continue
                
                if not reaction.products:
                    print("Error: No products found. Please check the equation format.")
                    continue
                
                reactions.append(reaction)
                
                print(f"Added reaction: {reaction.equation}")
                print(f"   Reactants: {dict(reaction.reactants)}")
                print(f"   Products: {dict(reaction.products)}")
                print(f"   Time: {reaction.time} minutes")
                break  
                
            except Exception as e:
                print(f"Error creating reaction: {e}")
                print("Please try again with correct format")
                continue

def input_compounds(reactions: List[Reaction]) -> Dict[str, float]:
    print("\n" + "="*70 + "\n")
    print("ENTER AVAILABLE COMPOUNDS")
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
            if compound_moles:
                break
            else:
                print("You must enter at least one compound before finishing.")
                continue
        
        if not compound:
            print("Please enter a valid compound name")
            continue
        
        amount = None
        while amount is None:
            amount_input = input(f"Amount of {compound}: ").strip()
            if not amount_input:
                print("Please enter an amount")
                continue
            try:
                amount = float(amount_input)
                if amount < 0:
                    print("Amount must be non-negative")
                    amount = None
                    continue
            except ValueError:
                print("Please enter a valid number")
                continue
        
        unit = None
        while unit is None:
            unit_input = input("Unit (g, kg, mol): ").strip().lower()
            if unit_input in ['g', 'gram', 'grams', 'kg', 'kilogram', 'kilograms', 'mol', 'mole', 'moles']:
                unit = unit_input
            else:
                print("Please enter a valid unit: g, kg, or mol")
                continue
        
        try:
            moles = convert_to_moles(amount, unit, compound)
            compound_moles[compound] = moles
            print(f"Added {moles:.3f} mol of {compound}")
        except Exception as e:
            print(f"Error adding compound: {e}")
            print("Please try again")
    
    print(f"\nAvailable compounds (in moles):")
    for compound, moles in sorted(compound_moles.items()):
        print(f"   {compound}: {moles:.3f} mol")
    
    return compound_moles

def input_constraints(reactions: List[Reaction]) -> Tuple[float, str]:
    print("\n" + "="*70 + "\n")
    print("ENTER OPTIMIZATION PARAMETERS")
    
    all_products = set()
    for reaction in reactions:
        all_products.update(reaction.products.keys())
    
    if all_products:
        print(f"Available products from your reactions: {', '.join(sorted(all_products))}")
    
    target = None
    while target is None:
        target_input = input("Enter target product to maximize: ").strip()
        if target_input:
            target = target_input
        else:
            print("Please enter a target product!")
    
    max_time = None
    while max_time is None:
        try:
            max_time_input = input("Maximum time limit (minutes): ").strip()
            if not max_time_input:
                print("Please enter a time limit")
                continue
            max_time = float(max_time_input)
            if max_time <= 0:
                print("Time must be positive!")
                max_time = None
                continue
        except ValueError:
            print("Please enter a valid number!")
            continue
    
    print(f"\nTarget product: {target}")
    print(f"Time limit: {max_time} minutes")
    
    return max_time, target

def print_results(max_yield: float, selected_reactions: List[str], stats: Dict, 
                 max_time: float, target_product: str):
    print("\n" + "="*70)
    print("                         OPTIMIZATION RESULTS                     ")
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
        
        if target_product in MOLAR_MASS:
            mass_yield = max_yield * MOLAR_MASS[target_product]
            print(f"Yield in grams: {mass_yield:.2f} g")
        
    else:
        print("No feasible solution found within the given constraints.")
        print("Try increasing the time limit or check compound availability.")

def main():
    print("="*70)
    print("                     CHEMICAL REACTION YIELD OPTIMIZER                     ")
    print("                     USING BRANCH AND BOUND ALGORITHM                     ")
    print("="*70)
    
    print("\nChoose input method:")
    print("1. CLI input")
    print("2. File input")
    
    choice = input("Enter choice (1-2): ").strip()

    while choice not in ["1", "2"]:
        print("Invalid choice. Please enter 1 or 2.")
        choice = input("\nEnter choice (1-2): ").strip()
    
    if choice == "2":
        example = input("\nDo you need an example input file? (y/n): ").strip().lower()
        while example not in ['y', 'n', 'yes', 'no']:
            print("Invalid choice. Please enter 'y' or 'n'.")
            example = input("\nDo you need an example input file? (y/n): ").strip().lower()

        if example in ['y', 'yes']:
            generate_input_example()
            print("Example input file generated at '../test/example.txt'")

        filename = input("\nEnter input filename (e.g., input.txt): ").strip()
        reactions, compound_moles, max_time, target_product = load_from_file(filename)
        
        while not reactions:
            print("Failed to load file (file name or format invalid).")
            filename = input("\nEnter input filename (e.g., input.txt): ").strip()
            reactions, compound_moles, max_time, target_product = load_from_file(filename)
        
        initial_data = {'reactions': reactions, 'compounds': compound_moles}
    else:
        reactions = input_reactions()
        compound_moles = input_compounds(reactions)
        max_time, target_product = input_constraints(reactions)
        initial_data = {'reactions': reactions, 'compounds': compound_moles}
    
    max_yield, selected_reactions, stats = branch_and_bound(
        reactions, max_time, compound_moles, target_product
    )
    
    print_results(max_yield, selected_reactions, stats, max_time, target_product)
    
    save_choice = input("\nDo you want to save the solution to a file? (y/n): ").strip().lower()
    while save_choice not in ['y', 'n', 'yes', 'no']:
        print("Invalid choice. Please enter 'y' or 'n'.") 
        save_choice = input("\nDo you want to save the solution to a file? (y/n): ").strip().lower()

    if save_choice in ['y', 'yes']:
        output_filename = input("Enter output filename (e.g., solution.txt): ").strip()
        if output_filename == "" or output_filename == ".txt" or output_filename == ".":
            output_filename = "solution.txt"

        if not output_filename.endswith('.txt'):
            output_filename += '.txt'

        save_solution_to_file(output_filename, max_yield, selected_reactions, 
                             stats, max_time, target_product, initial_data)

if __name__ == "__main__":
    main()