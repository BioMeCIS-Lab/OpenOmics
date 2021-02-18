COUNT = "_count"


def make_trie(words):
    """
    Args:
        words:
    """
    root = {}
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})

        if COUNT not in current_dict:
            current_dict[COUNT] = 1
        else:
            current_dict[COUNT] += 1
    return root


def longest_common_prefix(strs):
    """
    Args:
        strs:
    """
    def traverse_trie(dictionary, prefix):
        if len(prefix) > 100:
            return prefix

        if len(dictionary.keys()) == 1 and COUNT not in dictionary.keys():
            letter = list(dictionary.keys())[0]
            return traverse_trie(dictionary[letter], prefix + letter)
        else:
            return prefix

    trie = make_trie(strs)
    lcp_branches = []
    for branch in trie:
        branch_lcp = traverse_trie(trie[branch], branch)
        lcp_branches.append(branch_lcp)

    return lcp_branches
