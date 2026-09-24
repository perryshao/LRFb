function [trajectory, id] = readIPtxt(ip_txt)
% READIPTXT  Read an IP trajectory text file and its filename-encoded IDs.
%   Writes a sibling *_v1 file after replacing decimal commas with periods.
%   Returns head, left hand, right hand and wrist coordinates in columns.
%   id = [person, gesture, session, recording].

fidr = fopen(ip_txt, 'rt');
fidw = fopen([ip_txt '_v1'], 'wt');
while ~feof(fidr)
    s = fgetl(fidr);
    s = strrep(s, ',', '.');
    fprintf(fidw, '%s\n', s);
end
fclose(fidr);
fclose(fidw);

file = sprintf([ip_txt '_v1']);
fp = fopen(file);
frame_count = fscanf(fp, '%d', 1);
column_num = fscanf(fp, '%d', 1);
A = zeros(frame_count, column_num);
for i = 1:frame_count
    A(i, :) = fscanf(fp, '%f', column_num);
end
fclose(fp);
Head = A(:, 1:3);
Hand_Left = A(:, 4:6);
Hand_Right = A(:, 7:9);
Wrist = A(:, 10:12);
personid = str2double(ip_txt(end - 12:end - 11));
gestureid = str2double(ip_txt(end - 9:end - 8));
sessionid = str2double(ip_txt(end - 6));
recordid = str2double(ip_txt(end - 4));
id = [personid gestureid sessionid recordid];
trajectory = [Head Hand_Left Hand_Right Wrist];
