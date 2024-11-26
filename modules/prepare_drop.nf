// What exactly is this part performing?

// HGNC ID map

process prepare_drop {
    input:
    path input_file

    output:
    path 'hello_out.txt'

    script:
    """
    prepare_drop.sh ${input_file} > "hello_out.txt"
    """

    stub:
    """
    touch hello_out.txt
    """
}
