import time
from io_handler import *
from optimizer import branch_and_bound

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
    
    start_time = time.perf_counter()
    max_yield, selected_reactions, stats = branch_and_bound(
        reactions, max_time, compound_moles, target_product
    )
    end_time = time.perf_counter()
    
    processing_time = (end_time - start_time) * 1000
    stats['processing_time'] = processing_time
    
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