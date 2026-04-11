def search_line(keyword, script, multi="", annotationSymbol="#") -> list:
    n = 0
    if multi == "yes":
        arr = []
    for line in script:
        r1 = line.split(annotationSymbol, 1)[0]
        if len(r1) > 0:
            if r1.find(keyword) != -1:
                if multi != "yes":
                    return (line, n)
                else:
                    arr.append([line, n])
        n += 1
    return arr if multi == "yes" else ["0 0 0", -1]


def extra_section(n, script, break_arr=[], annotationSymbol="#"):
    outscript = []
    if len(break_arr) == 0:
        for line in script[n:]:
            if len(line.strip().split()) == 0:
                break
            outscript.append(line)
    else:
        for line in script[n:]:
            r1 = line.split(annotationSymbol, 1)[0]
            if len(r1) > 0:
                # if r1.strip("[").strip("]").strip() in break_arr:
                try:
                    if r1.split()[1].strip() in break_arr:
                        break
                except:  # noqa
                    pass
                outscript.append(line)
    return outscript


def extra_info(kws, script, break_arr, annotationSymbol="#"):
    data_dict = {}
    for kw in kws:
        n = search_line(kw, script, annotationSymbol=annotationSymbol)[1]
        if n == -1:
            data_dict[kw] = "NONE"
        else:
            n += 1
            section = extra_section(n, script, break_arr=break_arr, annotationSymbol=annotationSymbol)
            # section = extra_section(n, script, break_arr=[kw], annotationSymbol=annotationSymbol)
            if len(section) == 0:
                data_dict[kw] = "NONE"
            else:
                data_dict[kw] = []
                for line in section:
                    if len(line.split()) > 0:
                        if line.strip()[0] != "#":
                            data_dict[kw].append([s for s in line.strip().split()])
    return data_dict
