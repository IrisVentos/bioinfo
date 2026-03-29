## Dynamic Programming — the core idea

DP is fundamentally about one thing: don't solve the same subproblem twice. It applies whenever a problem has two properties:

Optimal substructure — the answer to a big problem is built from answers to smaller versions of the same problem

Overlapping subproblems — those smaller versions repeat across different paths of computation

The canonical mental model: a DAG of dependencies
Think of your problem as a directed acyclic graph where each node is a subproblem, and edges point from "things I need" to "things that need me". DP fills this graph in topological order — by the time you compute node x, everything it depends on is already done.

## The philosophical shift

### It remembers answers you've already computed so you don't compute them again.

That's it. Everything else is just details.
The mental shift DP asks for is subtle but important.
Naive thinking: "I need to explore all possibilities and find the best/count them all."
DP thinking: "I only need the answer to each subproblem, not the details of how it was reached. So I'll store just the answer and move on."
This is why DP can count billions of peptides without ever listing a single one — you track the number of ways rather than the ways themselves. The table compresses an exponentially large space into a polynomial amount of information.

### The two main ingredients

1. A recurrence (the dependence of n to subproblems n-1,n-2...)
2. A base case (where does the recursion bottom out, aka start)

