#ifndef SUBSET_H
#define SUBSET_H



#include <iostream>

struct subset_map
{
    int parent;
    int rank;
};
typedef subset_map Subset_map;

// A utility function to find set of an element i
// (uses path compression technique)
static int Find(std::map<int, Subset_map>& subsets, int i )
{
    // find root and make root as parent of i (path compression)
    if (subsets[i].parent != i)
        subsets[i].parent = Find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

// A function that does union of two sets of x and y by rank
// returns (lower root, upper root)
static std::pair<int,int> Union( std::map<int, Subset_map>& subsets, int x, int y )
{
    int xroot = Find( subsets, x );
    int yroot = Find( subsets, y );
    std::cout << "xroot = " << xroot << std::endl;
    std::cout << "yroot = " << yroot << std::endl;

    // Attach smaller rank tree under root of high rank tree
    // (Union by Rank)
    if ( subsets[xroot].rank < subsets[yroot].rank )
    {
        subsets[xroot].parent = yroot;
        return std::pair<int,int>(xroot,yroot);
    }
    else if (subsets[xroot].rank > subsets[yroot].rank)
    {
        subsets[yroot].parent = xroot;
        return std::pair<int,int>(yroot,xroot);
    }
    else
    {
        // If ranks are same, then make the larger name 0-simplex the root and increment
        // its rank by one
        std::cout << "ranks are the same" << std::endl;
        int big = std::max(xroot, yroot);
        int small = std::min(xroot, yroot);
        /*
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
        return std::pair<int,int>(yroot,xroot);*/
        subsets[small].parent = big;
        subsets[big].rank++;
        return std::pair<int,int>(small,big);
    }
}

#endif // SUBSET_H
