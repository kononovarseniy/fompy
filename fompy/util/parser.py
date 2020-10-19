import re


class Tokenizer:
    def __init__(self, rules):
        # Using re as lexical analyzer is described here https://stackoverflow.com/a/2359619
        self._re = re.compile('|'.join(f'(?P<{k}>{v})' for k, v in rules.items()), re.VERBOSE)

    def tokenize(self, text):
        pos = 0
        while pos < len(text):
            m = self._re.match(text, pos)
            if not m:
                raise SyntaxError(f'Tokenizer stopped at position {pos}')
            pos = m.end()
            token = m.lastgroup
            value = m.group(token)
            yield token, value


class TokenList(list):
    def __init__(self, tokens):
        super().__init__(tokens)
        self._pos = 0

    @property
    def position(self):
        return self._pos

    @position.setter
    def position(self, pos):
        self._pos = pos

    def get(self, offset=0):
        if self._pos + offset >= len(self):
            return None, None
        return self[self._pos + offset]

    @property
    def eof(self):
        return self._pos >= len(self)
