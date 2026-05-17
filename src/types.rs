use std::collections::HashMap;

pub type CoordinateMap = HashMap<(usize, usize), [usize; 6]>;

pub type RefArr = Vec<[usize; 6]>;

pub type RefMap = HashMap<String, (RefArr, usize)>;

pub type MapTide = HashMap<String, CoordinateMap>;

pub type RefLengths = HashMap<String, usize>;
