import sys


def InputParser(file_path: str) -> str:
    with open(file_path, 'r') as file:
        return file.readline().rstrip()


def RotateString(text: str) -> str:
    return text[-1]+text[:-1]


def BurrowsWheeler(text: str) -> str:
    text_pool = []
    
    rotated = text
    for _ in range(len(text)):
        rotated = RotateString(rotated)
        text_pool.append(rotated)
    
    text_pool.sort()

    return ''.join([x[-1] for x in text_pool])


if __name__ == '__main__':
    text = InputParser(sys.argv[1])
    print(BurrowsWheeler(text))