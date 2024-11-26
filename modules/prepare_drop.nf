// What exactly is this part performing?

// HGNC ID map

process prepare_drop {
    input:
    path input_file

    output:
    path 'hello_out.txt'

    script:
    """
    bash hello.sh ${input_file} > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}
