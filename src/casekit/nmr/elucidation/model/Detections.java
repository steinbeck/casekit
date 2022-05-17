package casekit.nmr.elucidation.model;

import casekit.nmr.model.DataSet;
import lombok.*;

import java.util.List;
import java.util.Map;
import java.util.Set;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
public class Detections {

    private Map<Integer, List<Integer>> detectedHybridizations;
    private Map<Integer, Map<String, Map<Integer, Set<Integer>>>> detectedConnectivities;
    private Map<Integer, Map<String, Map<Integer, Set<Integer>>>> forbiddenNeighbors;
    private Map<Integer, Map<String, Map<Integer, Set<Integer>>>> setNeighbors;
    private Map<Integer, Set<Integer>> fixedNeighbors;
    private List<DataSet> functionalGroups;
}
