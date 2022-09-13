function merge(v1, v2){
    if (v1 == "NA" || v1 == v2) {
        return v2
    } else if (v2 == "NA") {
        return v1
    } else {
        return "NA"
    }
}
