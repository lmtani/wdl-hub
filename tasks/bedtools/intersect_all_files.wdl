version 1.0


task IntersectAll {
    input {
        Array[File] regions
        String output_basename = "intersected.bed"
        Boolean stub = false
    }

    command <<<
    set -e

    bedtools --version > version.txt
    if [ ~{stub} = true ]; then
        touch ~{output_basename}
        exit 0
    fi

    # Define your list of BED files as a space-separated string
    bed_files="~{sep=' ' regions}"

    # Convert the string to an array
    read -r -a bed_array <<< "$bed_files"

    # Check if at least two BED files are provided
    if [ "${#bed_array[@]}" -lt 2 ]; then
        echo "Error: At least two BED files are required."
        exit 1
    fi

    # Use the first BED file as the initial intersection file
    cp "${bed_array[0]}" temp_intersection.bed

    # Loop through each BED file in the array and intersect
    for file in "${bed_array[@]:1}"
    do
        # Intersect with the current temporary file and update it
        bedtools intersect -a temp_intersection.bed -b "$file" > temp_intersection_update.bed
        mv temp_intersection_update.bed temp_intersection.bed
    done

    # Move the final intersection to a new file
    mv temp_intersection.bed ~{output_basename}

    >>>

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
    }

    output {
        File intersected = "~{output_basename}"
        File version = "version.txt"
    }
}
