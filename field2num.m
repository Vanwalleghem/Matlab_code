function a = field2num(s, f)
%FIELD2NUM put the values of a field into an array.
%*  A = FIELD2NUM(S, F) extract all the values of the field F of struct S
%   and put them in the numerical array A.
%
%*  See also: field2cell.

N = numel(s);
if N==0
    a = NaN;
    return
end
n = numel(s(1).(f));

a = NaN(N, n);
for i = 1:numel(s)
    a(i,:) = s(i).(f)(:);
end