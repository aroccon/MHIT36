function e_k = energy(k_mag,k_fit,e_k_fit)

% Returns the energy per unit |k| at k = k_mag (the physiscal wave number),
% right now the interpolation is linear in loglog-space, no 2 power law
% is assumed for the large scales (see Pope)
i = 1;
N = length(k_fit);

%e_k = ((k_mag - k_fit(1))/(k_fit(2) - k_fit(1)) * (e_k_fit(2) - e_k_fit(1)) + e_k_fit(1)) ;
e_k = exp((log(k_mag) - log(k_fit(1)))/(log(k_fit(2)) - log(k_fit(1))) * (log(e_k_fit(2)) - log(e_k_fit(1))) + log(e_k_fit(1)));
%e_k = e_k_fit(1)*(k_mag/k_fit(1))^2;
while (k_mag > k_fit(i))
    i = i+1;
    if i > N
        %e_k = ((k_mag - k_fit(N-1))/(k_fit(N) - k_fit(N-1)) * (e_k_fit(N) - e_k_fit(N-1)) + e_k_fit(N-1));
        e_k = exp((log(k_mag) - log(k_fit(N-1)))/(log(k_fit(N)) - log(k_fit(N-1))) * (log(e_k_fit(N)) - log(e_k_fit(N-1))) + log(e_k_fit(N-1)));
        break;
    end
	%e_k = ((k_mag - k_fit(i-1))/(k_fit(i) - k_fit(i-1)) * (e_k_fit(i) - e_k_fit(i-1)) + e_k_fit(i-1));
    e_k = exp((log(k_mag) - log(k_fit(i-1)))/(log(k_fit(i)) - log(k_fit(i-1))) * (log(e_k_fit(i)) - log(e_k_fit(i-1))) + log(e_k_fit(i-1)));
end

end




