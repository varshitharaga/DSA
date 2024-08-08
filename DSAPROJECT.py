from collections import deque


#Defines a class to represent a single nucleotide in a gene sequence.
class Nucleotide:
    def __init__(self, value):
        self.value = value
        self.next_gene = None

# Defines a class to represent a linked list of gene sequences.
class GeneLinkedList:
    def __init__(self):
        self.head = None

    def insert(self, gene_sequence):
        new_gene = Nucleotide(gene_sequence) #create a new nucleotide instance
        if not self.head:
            self.head = new_gene
        else:
            current = self.head #iterating to find last nucleotide
            while current.next_gene:
                current = current.next_gene
            current.next_gene = new_gene

#defining a binary tree node
class GeneTreeNode:
    def __init__(self, gene_sequence):
        self.value = gene_sequence
        self.left = None
        self.right = None
        self.linked_list = GeneLinkedList()

class GeneBinaryTree:
    def __init__(self):
        self.root = None

    def insert(self, gene_sequence):
        new_node = GeneTreeNode(gene_sequence)
        if not self.root:
            self.root = new_node
        else:
            self._insert_level_order(new_node)

    #helper function for insert
    #time complexity O(n)
    #worst case has to iterate through the entire tree to find the last position to insert the node (n is the number of nodes)
    def _insert_level_order(self, new_node):
        queue = deque()
        queue.append(self.root)

        while queue:
            current = queue.popleft()

            if not current.left:
                current.left = new_node
                return
            else:
                queue.append(current.left)

            if not current.right:
                current.right = new_node
                return
            else:
                queue.append(current.right)


    def print_level_order(self, root):
        if not root:
            return
        queue = deque()
        queue.append(root)
        while queue:
            node = queue.popleft()
            print(node.value, end=" ")
            if node.left:
                queue.append(node.left)
            if node.right:
                queue.append(node.right)


    #time complexity will be O(1)
    def compare_genes(self, gene1, gene2):
        if gene1 == gene2:
            return "Genes are identical."
        else:
            return "Genes are different."
        
    #time complexity is O(n+min(gene1,gene2))=O(n)
    def compare_sequences(self, gene1, gene2):
        node1 = self.search(gene1)
        node2 = self.search(gene2)
        if node1 is None or node2 is None:
            return

        differences = []
        min_length = min(len(gene1), len(gene2))  # Ensuring index is within the index of shorter sequence
        for i in range(min_length):
            if gene1[i] != gene2[i]:
                differences.append(i)

        if len(gene1) != len(gene2):
            differences.append(min_length)  # Add the index where the lengths differ

        print(f"Differences between {gene1} and {gene2}: {differences}")


    #worst case traverse the entire tree:O(n)
    def search(self, gene_sequence):
        if not self.root:
            return None

        queue = deque()
        queue.append(self.root)

        while queue:
            current = queue.popleft()
            if current.value == gene_sequence:
                return current

            if current.left:
                queue.append(current.left)

            if current.right:
                queue.append(current.right)

        return None
    
    #time complexity:O(n^2)
    def find_evolutionary_path(self, gene):
        node = self.search(gene)
        if node is None:
            return

        path = []
        while node:
            path.append(node.value) 
            node = self.find_parent(node.value)

        print(f"Evolutionary path for {gene}: {' -> '.join(path[::-1])}")

    #helper function
    def find_parent(self, gene):
        return self._find_parent(gene, self.root)  


    #helper function
    def _find_parent(self, gene, node):
        if node is None or node.value == gene:
            return None

        if node.left and node.left.value == gene:
            return node
        elif node.right and node.right.value == gene:
            return node

        # Recursive call should prioritize going up the path containing the target gene
        parent = self._find_parent(gene, node.left)
        if parent:
            return parent
        return self._find_parent(gene, node.right)

#time complexity:O(n+m)
    def evolutionary_rate(self, gene1, gene2):
        node1 = self.search(gene1)
        node2 = self.search(gene2)
        if node1 is None or node2 is None:
            return None

        seq1 = node1.value
        seq2 = node2.value

        if len(seq1) != len(seq2):
            return None  # Sequences must be of equal length for rate calculation

        substitutions = sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])
        rate = substitutions / len(seq1)
        return rate

  #time complexity:O(n)
    def gene_relation(self, gene1, gene2):
        # Search for both genes in the binary tree
        node1 = self.search(gene1)
        node2 = self.search(gene2)

        # If either gene is not found, return None
        if node1 is None or node2 is None:
            return None

        # Find the common ancestor of both genes
        common_ancestor = self.find_common_ancestor(self.root, gene1, gene2)

        # Determine the relationship between gene1 and gene2 based on their positions relative to the common ancestor
        if common_ancestor.value == gene1 or common_ancestor.value == gene2:
            return f"{gene1} is an ancestor of {gene2}" if common_ancestor.value == gene1 else f"{gene2} is an ancestor of {gene1}"
        
    
    def find_common_ancestor(self, node, gene1, gene2):
        # Helper function to find the common ancestor of two genes
        if node is None or node.value == gene1 or node.value == gene2:
            return node

        left_ancestor = self.find_common_ancestor(node.left, gene1, gene2)
        right_ancestor = self.find_common_ancestor(node.right, gene1, gene2)

        if left_ancestor and right_ancestor:
            return node
        elif left_ancestor:
            return left_ancestor
        else:
            return right_ancestor

    #time complexity:O(n)
    def find_generation(self, gene):
        # Find the node corresponding to the gene
        node = self.search(gene)

        # If gene is not found, return None
        if node is None:
            return None

        # Find the generation of the gene
        generation = self._find_generation(self.root, gene, 1)
        return generation

    def _find_generation(self, node, gene, current_generation):
        if node is None:
            return None  # Gene not found

        if node.value == gene:
            return current_generation

        left_generation = self._find_generation(node.left, gene, current_generation + 1)
        if left_generation is not None:
            return left_generation

        right_generation = self._find_generation(node.right, gene, current_generation + 1)
        return right_generation



# Example usage
# ATCGCTAGCTAG
# ATCGATTGATCG
# ATCGCGTACGTA
# ATCGCGTACGAT
# ATCGCGTACTGA
# ATCGCGTACAGT
# ATCGCGTACATG

if __name__ == "__main__":
    binary_tree = GeneBinaryTree()

    # Ask the user for the number of gene sequences to insert
    num_sequences = int(input("Enter the number of gene sequences to insert: "))
    for _ in range(num_sequences):
        gene_sequence = input("Enter the gene sequence to insert: ")
        binary_tree.insert(gene_sequence)

    # Automatically print the level order traversal
    print("Level order traversal:")
    binary_tree.print_level_order(binary_tree.root)
    print()

    while True:
        print("\nMenu:")
        print("1. Compare genes")
        print("2. Find evolutionary path")
        print("3. Compare sequences")
        print("4. Calculate evolutionary rate")
        print("5. Determine gene relation")
        print("6. Find generation")
        print("7. Exit")

        choice = input("\nEnter your choice: ")

        if choice == "1":
            gene1 = input("Enter the first gene sequence: ")
            gene2 = input("Enter the second gene sequence: ")
            print(binary_tree.compare_genes(gene1, gene2))

        elif choice == "2":
            gene = input("Enter the gene sequence: ")
            binary_tree.find_evolutionary_path(gene)

        elif choice == "3":
            gene1 = input("Enter the first gene sequence: ")
            gene2 = input("Enter the second gene sequence: ")
            binary_tree.compare_sequences(gene1, gene2)

        elif choice == "4":
            gene1 = input("Enter the first gene sequence: ")
            gene2 = input("Enter the second gene sequence: ")
            print("Evolutionary rate:", binary_tree.evolutionary_rate(gene1, gene2))

        elif choice == "5":
            gene1 = input("Enter the first gene sequence: ")
            gene2 = input("Enter the second gene sequence: ")
            print(binary_tree.gene_relation(gene1, gene2))

        elif choice == "6":
            gene = input("Enter the gene sequence: ")
            print("Generation:", binary_tree.find_generation(gene))

        elif choice == "7":
            print("Exiting...")
            break

        else:
            print("Invalid choice. Please try again.")
