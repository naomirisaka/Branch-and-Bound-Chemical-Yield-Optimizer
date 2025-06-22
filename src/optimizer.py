import heapq
from typing import Dict, List, Tuple
from reactions import Reaction, Node

def can_apply(reaction: Reaction, available_moles: Dict[str, float]) -> bool:
    return all(reactant in available_moles and available_moles[reactant] >= coeff
               for reactant, coeff in reaction.reactants.items())

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
                projected_time = time + actual_time
                if projected_time <= max_time:
                    time = projected_time
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
        return potential_yield / actual_time if actual_time > 0 else 0
    
    reactions.sort(key=efficiency, reverse=True)

    queue = []
    root = Node(-1, 0, 0, initial_moles, 0, [])
    root.bound = calculate_bound(root, reactions, max_time, target_product)
    heapq.heappush(queue, root)

    max_yield = 0
    best_combination = []
    nodes_explored = 0

    while queue:
        node = heapq.heappop(queue)
        nodes_explored += 1
        
        # prune if bound is not better than current best
        if node.bound <= max_yield:
            continue
            
        if node.level == len(reactions) - 1:
            continue

        i = node.level + 1
        reaction = reactions[i]

        # xi = 1, take the current reaction
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
                     
                    if taken_node.bound > max_yield:
                        heapq.heappush(queue, taken_node)

        # xi = 0, skip the current reaction
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

    final_yield = current_moles.get(target_product, 0)
    stats = {
        'nodes_explored': nodes_explored,
        'total_time': total_time,
        'efficiency': final_yield / max(total_time, 0.1)
    }
    return final_yield, selected_reactions, stats