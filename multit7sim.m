n = [470,470+500,470+1000,470+1500,470+2000];
for i = 1:length(n)
    coupled_ode_run(n(i))
end
