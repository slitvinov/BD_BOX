function [delta_array, array_corr] = polymer_b2bcorr(R)
% delta_array: array of differences between beads
% array_corr: array of correlation coefficients

number_of_beads = size(R, 3);

% start from this bead of the polymer
first_bead=4;
last_bead=number_of_beads-first_bead;

% delta
delta_array = 0:(last_bead-first_bead);
array_corr = zeros(size(delta_array));
for dn=1:length(delta_array)
  di = delta_array(dn);
  % initial bead
  sum_corr = 0;
  nsum = 0;
  for initial_bead=first_bead:(last_bead-di)
     next_bead = initial_bead + di;
     for dim=1:3
       coord1 = R(:, dim, initial_bead);
       coord2 = R(:, dim, next_bead);
       sum_corr = sum_corr + corrcoef(coord1, coord2);
       nsum = nsum + 1;
     end
   end
   array_corr(dn) = sum_corr/nsum;
end
endfunction