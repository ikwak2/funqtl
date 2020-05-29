# determine cross type
crosstype <-
    function(cross)
    {
        type <- class(cross)
        type <- type[type != "cross" & type != "list"]
        if(length(type) > 1) {
            warning("cross has multiple classes")
        }
        type[1]
    }
