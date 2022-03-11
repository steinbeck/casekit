package casekit.nmr.model.nmrium;

import casekit.nmr.model.PathLength;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import lombok.*;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
@JsonIgnoreProperties(ignoreUnknown = true)
public class J {

    private PathLength pathLength;
}
