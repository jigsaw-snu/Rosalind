import sys


'''
    <InputParser>
    Parse Input file and save data into dictionary
'''
def InputParser(input_path: str) -> dict:
    params = dict()

    with open(input_path, 'r') as file:
        # first line implies for total amount of money
        line = file.readline().rstrip()
        params['money'] = int(line)

        # second line implies for list of available coins
        line = file.readline().rstrip()

        if ',' in line:
            params['coins'] = list(map(int, line.split(',')))
        else:
            params['coins'] = list(map(int, line.split()))

    return params


'''
    <MinNumCoins>
    Find Minimum Number of Coins with given money
'''
# dictionary for memoization
ledger = {0:0}

def MinNumCoins(money: int, coins: list) -> int:
    if not money:
        return ledger[0]

    for m in range(1, money+1):
        num_coins = [ledger[m-c] for c in coins if c <= m]
        ledger[m] = min(num_coins)+1

    return ledger[money]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(MinNumCoins(params['money'], params['coins']))
