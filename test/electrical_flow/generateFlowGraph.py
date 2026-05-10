import random
import argparse

def generateFlow(N, M, U):
    all_edges = [(u, v) for u in range(N) for v in range(u+1, N)]
    edges = random.sample(all_edges, M)
    weighted_edges = [(u, v, random.randint(1, U)) for u, v in edges]
    random.shuffle(weighted_edges)
    return weighted_edges


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='generatePerfectGraph',
        description='simple graph generator for graphs with a perfect matching'
    )
    parser.add_argument("N", type=int, help='Number of vertices')
    parser.add_argument("M", type=int, help='Number of edges')
    parser.add_argument("U", type=int, help='Max capacity of an edge')
    args = parser.parse_args()

    edges = generateFlow(args.N, args.M, args.U)

    print('p', 'max', args.N, len(edges))
    print('n 1 s')
    print(f'n {args.N} t')
    for (u, v, c) in edges:
        print('e', u+1, v+1, c)

