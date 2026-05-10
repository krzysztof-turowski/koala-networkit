import random
import argparse



def generateGeneralGraph(N, M, isBipartite):
    assert N%2 == 0, 'N must be even'
    assert M >= N//2, 'M must not be smaller than N/2'

    vertices = [v for v in range(N)]
    random.shuffle(vertices)

    perfect = {(min(vertices[i:i+2]), max(vertices[i:i+2])) for i in range(0, N, 2)}

    all_edges = {(u, v) for u in range(N) for v in range(u+1, N)}
    all_edges.difference(perfect)

    edges = list(perfect)
    edges += random.sample(all_edges, M-len(edges))
    random.shuffle(edges)
    
    return edges


def generateBipartiteGraph(N, M, isBipartite):
    assert N%2 == 0, 'N must be even'
    assert M >= N//2, 'M must not be smaller than N/2'

    V = random.sample([v for v in range(N)], N//2)
    U = list({v for v in range(N)}.difference(V))
    random.shuffle(U), random.shuffle(V)

    perfect = set(zip(U, V))

    all_edges = {(u, v) for u in U for v in V}
    all_edges.difference(perfect)

    edges = list(perfect)
    edges += random.sample(all_edges, M-len(edges))
    random.shuffle(edges)

    return edges


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='generatePerfectGraph',
        description='simple graph generator for graphs with a perfect matching'
    )
    parser.add_argument("N", type=int, help='Number of vertices')
    parser.add_argument("M", type=int, help='Number of edges')
    parser.add_argument('-bp', '--bipartite', type=bool,
                        default=False, help='Should graph be bipartite')
    args = parser.parse_args()

    if args.bipartite:
        edges = generateBipartiteGraph(args.N, args.M, args.bipartite)
    else:
        edges = generateGeneralGraph(args.N, args.M, args.bipartite)

    print('p', 'edge', args.N, len(edges))
    for (u, v) in edges:
        print('e', u+1, v+1)
