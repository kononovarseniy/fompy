"""
This module contains classes for lexical analysis of text.
"""

import re


class Tokenizer:
    """
    A class to parse text.
    """

    def __init__(self, rules):
        """
        Initialize the parser.

        Parameters
        ----------
        rules : dict(str 'token_type' : str 'regex')
            A set of parsing rules establishing the correspondence between token types and regular expressions.
        """
        # Using re as lexical analyzer is described here https://stackoverflow.com/a/2359619
        self._re = re.compile('|'.join(f'(?P<{k}>{v})' for k, v in rules.items()), re.VERBOSE)

    def tokenize(self, text):
        """
        Parse `text` and generate tokens.

        Parameters
        ----------
        text : str

        Returns
        -------
        str 'token', str 'value'
            The current token.
        Raises
        ------
        SyntaxError
            The parser has failed to recognize a token.
        """
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
    """
    A class functioning as a list of tokens (extends `list`).
    """

    def __init__(self, tokens):
        super().__init__(tokens)
        self._pos = 0

    @property
    def position(self):
        """Get or set the current position."""
        return self._pos

    @position.setter
    def position(self, pos):
        self._pos = pos

    def get(self, offset=0):
        """
        Get a token.

        Parameters
        ----------
        offset : int
            An offset from the current `position`.

        Returns
        -------
        `self[ position + offset ]`
        None, None
            When `position + offset` is past the end of the token list.
        """
        if self._pos + offset >= len(self):
            return None, None
        return self[self._pos + offset]

    @property
    def eof(self):
        """Get a boolean that is `True` if `position` is past the end of the token list, `False` otherwise."""
        return self._pos >= len(self)
