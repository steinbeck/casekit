package casekit.nmr.model.nmrium;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import lombok.*;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
@JsonIgnoreProperties(ignoreUnknown = true)
public class Signal {

    private String id;
    private String kind;
    private String multiplicity;
    private Integer sign;
}
