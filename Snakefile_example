parameter = [1, 2, 3]

rule make_firstfile:
    output:
        'firstfile_{filename}.txt'
    shell:
        'uname -a > firstfile_{wildcards.filename}.txt'

rule make_secondfile:
    input:
        'firstfile_{filename}.txt'
    output:
        'secondfile_{filename}.txt'
    shell:
        'cat firstfile_{wildcards.filename}.txt > secondfile_{wildcards.filename}.txt'

