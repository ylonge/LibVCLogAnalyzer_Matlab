function [qstep] = qp2qstep(qp)
qstep = 2.^((qp - 4)/6);
end