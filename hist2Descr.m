function [feat] = hist2Descr(feat,descr,descr_mag_thr)
% Function: 将最后一步中的直方图转换为描述符（特征向量）
descr = descr/norm(descr);
descr = min(descr_mag_thr,descr);
descr = descr/norm(descr);
feat.descr = descr;
end