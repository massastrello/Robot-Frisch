function IDX = cost_int_time7(x)
    IDX = abs(x(1)) + norm(x(2:3));%+ norm(x(2:3)); 
end