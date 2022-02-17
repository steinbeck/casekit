package casekit.nmr.elucidation.model;

import lombok.*;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

@AllArgsConstructor
@NoArgsConstructor
@Getter
@Setter
@ToString
public class Grouping {

    Map<String, Double> tolerances = new HashMap<>();
    Map<String, Map<Integer, List<Integer>>> groups;
    Map<String, Map<Integer, Integer>> transformedGroups;
}
