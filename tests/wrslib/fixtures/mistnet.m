function [preds, probs, classes, y, x, elevs] = mistnet(~)
%MISTNET Deterministic Cartesian prediction grid for remapping tests.

y = [-1000 0 1000];
x = [-1000 0 1000];
elevs = [-1 1];
base_predictions = [11 12 13; 21 22 23; 31 32 33];
preds = cat(3, base_predictions, base_predictions);
probs = cat(4, preds, preds + 100);
classes.rain = 3;
end
