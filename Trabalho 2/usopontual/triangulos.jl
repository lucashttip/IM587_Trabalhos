function racos(a, b, c)
    c = -(a^2 - b^2 - c^2)/(2*b*c)
    return acosd(c)
end