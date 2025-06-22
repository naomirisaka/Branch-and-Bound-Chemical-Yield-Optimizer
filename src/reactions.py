import re
from typing import Dict, List, Tuple

# periodic table data
ATOMIC_MASS = {
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

MOLAR_MASS = {}

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
            
            # parse left side (reactants)
            for r in lhs.strip().split('+'):
                coeff, compound = self.extract_compound_with_coefficient(r.strip())
                if compound:
                    reactants[compound] = coeff
            
            # parse right side (products)
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
                return 0.0 
        
        return max_extent if max_extent != float('inf') else 0.0

    def __str__(self):
        return f"{self.equation} | time={self.time} min"

class Node:
    def __init__(self, level: int, total_yield: float, total_time: float, 
                 available_moles: Dict[str, float], bound: float, taken: List[bool]):
        self.level = level
        self.total_yield = total_yield
        self.total_time = total_time
        self.available_moles = available_moles.copy()
        self.bound = bound
        self.taken = taken.copy()

    def __lt__(self, other):
        return self.bound > other.bound

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
    if formula in MOLAR_MASS:
        return MOLAR_MASS[formula]
    
    try:
        element_counts = parse_chemical_formula(formula)
        total_weight = 0.0
        
        for element, count in element_counts.items():
            if element in ATOMIC_MASS:
                total_weight += ATOMIC_MASS[element] * count
            else:
                print(f"Unknown element: {element}")
                while True:
                    try:
                        weight = float(input(f"Enter atomic mass for {element} (g/mol): "))
                        if weight > 0:
                            ATOMIC_MASS[element] = weight
                            total_weight += weight * count
                            break
                        print("Atomic mass must be positive!")
                    except ValueError:
                        print("Please enter a valid number!")
        
        MOLAR_MASS[formula] = total_weight 
        return total_weight
        
    except Exception as e:
        print(f"Error calculating molecular weight for {formula}: {e}")
        while True:
            try:
                mw = float(input(f"Enter molecular weight for {formula} (g/mol): "))
                if mw > 0:
                    MOLAR_MASS[formula] = mw
                    return mw
                print("Molecular weight must be positive.")
            except ValueError:
                print("Please enter a valid number.")

def convert_to_moles(amount: float, unit: str, compound: str) -> float:
    unit = unit.lower().strip()
    
    if unit in ['mol', 'moles']:
        return amount
    elif unit in ['g', 'gram', 'grams']:
        mw = calculate_molecular_weight(compound)
        return amount / mw
    elif unit in ['kg', 'kilogram', 'kilograms']:
        mw = calculate_molecular_weight(compound)
        return (amount * 1000) / mw
    else:
        print(f"Unknown unit: {unit}")
        print("Supported units: mol, g, kg")
        print("Assuming grams...")
        mw = calculate_molecular_weight(compound)
        return amount / mw